#' Create a nnspline object
#'
#' @param x A numeric matrix giving the locations (as rows) for which the spline values are desired
#' @param nodes A numeric matrix giving the spline knot locations as rows.
#' @param n_parents The number of parents each node should have.
#' @param parameters The parameters used in the covariance function.
#' @param covariance_function A function used to compute covariances between spline inputs. Should have the arguments x1 and x2, vectors giving the spline inputs, and p, a vector of parameters.
#' @param node_graph An optional argument. If supplied, it needs to be an igraph object representing a directed acyclic graph.
#'
#' @return A multispline object as a list with the following elements:
#'   - values The output values of the spline at x.
#'   - x The locations at which the spline is evaluated.
#'   - x_parents The nodes which are parents of each x.
#'   - node_values The output values of the spline at the nodes.
#'   - node_values_map Mapping indices for the node values. See ?TMB::MakeADFun.
#'   - nodes The locations of the spline nodes.
#'   - node_graph An igraph directed acyclic graph represnting the node parent structure.
#'   - parameters The parameters for the covariance function.
#'   - parameter_map Mapping indices for the parameters. See ?TMB::MakeADFun.
#'   - covariance_function The function used to compute covariances between input points.
#'   - x_covariance Sparse covariance matrices for the spline x values.
#'   - node_covariance Sparse covariance matrices for the spline nodes.
#' 
#' @export
create_nnspline<- function(
    x,
    nodes = x,
    n_parents = 4,
    parameters = c(1, 1),
    covariance_function = function(x1, x2, p) {
      d<- sqrt(sum((x1 - x2)^2))
      variance<- p[[1]]
      range<- p[[2]]
      cov<- variance * (1 + d / range) * exp( -d / range )
      return(cov)
    },
    node_graph
  ) {
  if( !("matrix" %in% class(x)) ) x<- matrix(x, ncol = 1)
  if( !("matrix" %in% class(nodes)) ) nodes<- matrix(nodes, ncol = 1)

  t_nodes<- t(nodes)
  x_parents<- t(
    apply(
        x,
        MARGIN = 1,
        function(row) {
          d<- t_nodes - row
          d<- sqrt(colSums(d^2))
          if( any( d < .Machine$double.eps ) ) {
            p<- which( d < .Machine$double.eps )[[1]] # There shouldn't be more than one, but just to be safe
            d[-p]<- 0
            d[p]<- -1
          } else {
            p<- order(d)[1:n_parents]
            p<- p[!is.na(p)]
            d[-p]<- 0
            d[p]<- 1
          }
          return(d)
        }
    )
  )

  x_parents<- Matrix::Matrix(x_parents, sparse = TRUE, doDiag = FALSE)
  x_covariance<- x_parents

  ind<- (x_parents@x == 1)
  mat_i<- x_covariance@i[ind]
  mat_i<- c(mat_i, unique(mat_i))
  mat_j<- c(
    ip_to_j(x_covariance@i, x_covariance@p)[ind] + 1,
    rep(0, length(unique(mat_i)))
  )
  mat_x<- c(
    x_covariance@x[ind],
    rep(1, length(unique(mat_i)))
  )
  x_covariance<- Matrix::sparseMatrix(
    i = mat_i + 1,
    j = mat_j + 1,
    x = mat_x
  )

  if( missing(node_graph) ) {
    node_graph<- distance_matrix_to_dag(
        dist(nodes),
        k = n_parents
    )
  }
  node_pairs<- get_pairs(node_graph, x_parents)
  node_covariance<- Matrix::sparseMatrix(
      i = node_pairs$i,
      j = node_pairs$j,
      x = rep(1, nrow(node_pairs))
    )

  spline<- list(
    values = numeric(nrow(x)),
    x = x,
    x_parents = x_parents,
    node_values = numeric(nrow(nodes)),
    node_values_map = seq(nrow(nodes)),
    nodes = nodes,
    node_graph = node_graph,
    parameters = parameters,
    parameter_map = seq_along(parameters),
    covariance_function = covariance_function,
    x_covariance = x_covariance,
    node_covariance = node_covariance
  )
  spline<- update_spline_covariance(
    spline,
    spline$parameters
  )
  spline<- update_spline_values(
    spline,
    spline$node_values
  )
  return(spline)
}

#' @describeIn update_spline Update a spline's covariance matrix
#' 
#' @param parameters The new parameters for the spline. If missing, will
#'   use the pre-existing parameters stored in the spline.
#' @param only_node_covariance Logical. If TRUE, the cross-covariance matrix between x and nodes will not be updated.
#' 
#' @export
update_spline_covariance<- function(
    spline,
    parameters,
    only_node_covariance = FALSE,
    ...
  ) {
  if( !missing(parameters) ) spline$parameters<- parameters
  mode<- "numeric"
  if( "advector" %in% class(spline$parameters) && requireNamespace("RTMB") ) mode<- "advector"
  if( mode == "advector" ) {
    RTMB::ADoverload()
    spline$x_covariance<- as(spline$x_covariance, "adsparse")
    spline$node_covariance<- as(spline$node_covariance, "adsparse")
  } else {}

  node_m<- cbind(
    spline$node_covariance@i + 1,
    ip_to_j(spline$node_covariance@i, spline$node_covariance@p) + 1
  )
  covs<- apply(
    node_m,
    MARGIN = 1,
    function(row) {
      cov<- spline$covariance_function(
        spline$nodes[row[[1]], ],
        spline$nodes[row[[2]], ],
        spline$parameters
      )
      return( cov )
    },
    simplify = FALSE
  )
  spline$node_covariance@x<- do.call(c, covs)

  if( only_node_covariance ) return(spline)

  if( any(spline$x_parents@x != -1) ) {
    x_m<- cbind(
      spline$x_covariance@i + 1,
      ip_to_j(spline$x_covariance@i, spline$x_covariance@p) + 1
    )
    covs<- apply(
      x_m,
      MARGIN = 1,
      function(row) {
        if( row[[2]] == 1 ) {
          loc2<- spline$x[row[[1]], ]
        } else {
          loc2<- spline$nodes[row[[2]] - 1, ]
        }
        cov<- spline$covariance_function(
          spline$x[row[[1]], ],
          loc2,
          spline$parameters
        )
        return(cov)
      },
      simplify = FALSE
    )
    spline$x_covariance@x<- do.call(c, covs)
  } else {}

  return(spline)
}

#' @describeIn update_spline Update a spline's values
#' 
#' @param node_values The new new values for the spline. If missing, will use
#'   the pre-existing node values stored in the spline.
#' 
#' @export
update_spline_values<- function(
    spline,
    node_values,
    ...
  ) {
  if( !missing(node_values) ) spline$node_values<- node_values
  mode<- "numeric"
  if( "advector" %in% class(spline$node_values) && requireNamespace("RTMB") ) mode<- "advector"
  if( mode == "advector" ) spline$values<- RTMB::advector(spline$values)

  x_type<- apply(
    spline$x_parents,
    MARGIN = 1,
    function(parents) any(parents == -1)
  )
  # Update x values that are actually node values
  node_x<- which(x_type)
  node_x_parents<- apply(
    spline$x_parents[node_x, , drop = FALSE],
    MARGIN = 1,
    function(parents) which(parents == -1)
  )
  spline$values[node_x]<- spline$node_values[node_x_parents]

  # Update x values that are between nodes
  no_node_x<- which(!x_type)
  if( length(no_node_x) > 0 ) {
    no_node_parents<- apply(
      spline$x_parents[no_node_x, ],
      MARGIN = 1,
      function(parents) return(which(parents == 1)),
      simplify = FALSE
    )
    no_node_prediction<- lapply(
      seq_along(no_node_x),
      function(i) {
        p<- no_node_parents[[i]]
        Sigma<- spline$node_covariance[c(1, p), c(1, p)]
        if( mode == "advector" ) {
          Sigma<- adsparse_to_admatrix(Sigma)
        } else {
          Sigma<- as.matrix(Sigma)
        }
        Sigma[1, ]<- spline$x_covariance[no_node_x[i], c(1, p + 1)]
        Sigma[-1, 1]<- 0
        Sigma<- Sigma + t(Sigma)
        diag(Sigma)<- diag(Sigma) / 2

        cmvn<- conditional_normal(
          cbind(c(spline$values[i], spline$node_values[p])),
          cbind(numeric(length(p) + 1)),
          Sigma,
          predicted = 1,
          predictors = seq_along(p) + 1
        )
        return(cmvn$mean)
      }
    )
    spline$values[no_node_x]<- do.call(c, no_node_prediction)
  } else {}
  
  return(spline)
}


#' Update a spline based on node values and parameters values
#'
#' @param spline An nnspline
#' @param ... Used to pass arguments onto sub-functions.
#'
#' @return An updated copy of the spline.
#'
#' @export
update_spline<- function(
    spline,
    ...
  ) {
  spline<- update_spline_covariance(spline, ...)
  spline<- update_spline_values(spline, ...)

  return(spline)
}

#' Evaluate the log-likelihood of a spline
#'
#' @param x The node values
#' @param spline The spline object
#' @param log If TRUE, returns the log-likelihood.
#'
#' @return The log-likelihood of the spline evaluated at x. If x is a simref
#'   object, then the values of x will be replaced with simulated values.
#'
#' @export
dspline<- function(x, spline, log = TRUE) {
  mode<- "numeric"
  if( ("advector" %in% class(x) || "adsparse" %in% class(spline$node_covariance)) && requireNamespace("RTMB") ) mode<- "advector"
  node_order<- as.numeric(igraph::topo_sort(spline$node_graph))
  ll<- 0
  for( i in node_order ) {
    p<- as.numeric(igraph::neighbors(spline$node_graph, i, "in"))
    Sigma<- spline$node_covariance[c(i, p), c(i, p), drop = FALSE]
    if( mode == "advector" ) {
      Sigma<- adsparse_to_admatrix(Sigma)
    } else {
      Sigma<- as.matrix(Sigma)
    }
    Sigma<- Sigma + t(Sigma)
    diag(Sigma)<- diag(Sigma) / 2
    cmvn<- conditional_normal(
      cbind(x[c(i, p)]),
      cbind(0 * spline$node_values[c(i, p)]),
      Sigma
    )
    if( mode == "advector" | "simref" %in% class(x) ) {
      ll<- ll + RTMB::dmvnorm(
        x[i],
        cmvn$mean,
        cmvn$covariance,
        log = TRUE
      )
    } else {
      ll<- ll + mvtnorm::dmvnorm(
        x[i],
        cmvn$mean,
        cmvn$covariance,
        log = TRUE
      )
    }
  }
  if( !log ) ll<- exp(ll)
  return(ll)
}

#' Simulate a spline
#' 
#' @param n The number of samples
#' @param spline The spline to simulate from
#' 
#' @return A list of length n with the simulated splines
#' 
#' @export
rspline<- function(
    n = 1,
    spline
  ) {
  spline<- update_spline_covariance(spline)
  node_values<- matrix(
    0,
    nrow = length(spline$node_values),
    ncol = n
  )
  node_order<- as.numeric(igraph::topo_sort(spline$node_graph))
  for( i in node_order ) {
    p<- as.numeric(igraph::neighbors(spline$node_graph, i, "in"))
    Sigma<- as(
      spline$node_covariance[c(i, p), c(i, p)],
      "matrix"
    )
    Sigma<- Sigma + t(Sigma)
    diag(Sigma)<- diag(Sigma) / 2
    cmvn<- conditional_normal(
      node_values[c(i, p), , drop = FALSE],
      0 * node_values[c(i, p), , drop = FALSE],
      Sigma
    )
    node_values[i, ]<- MASS::mvrnorm(
      n = ncol(node_values),
      mu = 0 * cmvn$mean[, 1],
      Sigma = cmvn$covariance
    )
    node_values[i, ]<- node_values[i, ] + cmvn$mean
  }
  node_values<- sweep(
    node_values,
    MARGIN = 2,
    STATS = colMeans(node_values)
  )
  sims<- lapply(
    seq(n),
    function(i) update_spline_values(spline, node_values[, i])
  )

  return(sims)
}





#' Evaluate a nnspline
#' 
#' @param x The input values
#' @param spline The spline
#' @param index Logical. If TRUE, returns the index instead of the value.
#' 
#' @return The value(s) of the spline at x
#'
#' @export
nns<- function(
    x,
    spline,
    index = FALSE
  ) {
  if( !("matrix" %in% class(x)) ) x<- cbind(x)
  spline_x<- t(spline$x)
  idx<- apply(
    x,
    MARGIN = 1,
    function(row) {
        d<- spline_x - row 
        d<- sqrt(colSums(d^2))
        idx<- which.min(d)
        return(idx)
    }
  )
  if( index ) return(idx) else return(spline$values[idx])
}

