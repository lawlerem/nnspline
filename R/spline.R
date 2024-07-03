#' Create a nnspline object
#'
#' @param x A numeric matrix giving the locations (as rows) for which the spline values are desired
#' @param nodes A numeric matrix giving the spline knot locations as rows.
#' @param n_parents The number of parents each node should have.
#' @param variance The marginal variance of the spline.
#' @param parameters The parameters used in the correlation function.
#' @param correlation_function A function used to compute correlation between spline inputs. Should have the arguments x1 and x2, vectors giving the spline inputs, and p, a vector of parameters.
#'
#' @return A multispline object as a list with the following elements:
#'   - values The output values of the spline at x.
#'   - x The locations at which the spline is evaluated.
#'   - x_parents The nodes which are parents of each x.
#'   - node_values The output values of the spline at the nodes.
#'   - node_values_map Mapping indices for the node values. See ?TMB::MakeADFun.
#'   - nodes The locations of the spline nodes.
#'   - node_graph An igraph directed acyclic graph represnting the node parent structure.
#'   - variance The marginal variance of the spline.
#'   - parameters The parameters for the correlation function.
#'   - parameter_map Mapping indices for the parameters. See ?TMB::MakeADFun.
#'   - correlation_function The function used to compute correlations between input points.
#'   - x_covariance Sparse covariance matrices for the spline x values.
#'   - node_covariance Sparse covariance matrices for the spline nodes.
#' 
#' @export
create_nnspline<- function(
    x,
    nodes = x,
    n_parents = 4,
    variance = 1,
    parameters = c(1),
    correlation_function = function(x1, x2, p) {
      d<- sqrt(sum((x1 - x2)^2))
      range<- p[[1]]
      (1 + d / range) * exp( -d / range )
    }
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

  x_parents<- Matrix::Matrix(x_parents, sparse = TRUE)
  x_covariance<- x_parents

  ind<- (x_parents@x == 1)
  mat_i<- x_covariance@i[ind]
  mat_j<- nnspline:::ip_to_j(x_covariance@i, x_covariance@p)[ind]
  mat_x<- x_covariance@x[ind]
  x_covariance<- Matrix::sparseMatrix(
    i = mat_i + 1,
    j = mat_j + 1,
    x = mat_x
  )

  node_graph<- nnspline:::distance_matrix_to_dag(
    dist(nodes),
    k = n_parents
  )
  node_pairs<- nnspline:::get_pairs(node_graph, x_parents)
  node_covariance<- Matrix::sparseMatrix(
      i = node_pairs$i,
      j = node_pairs$j,
      x = rep(1, nrow(node_pairs))
    )

  return(
    list(
      values = numeric(nrow(x)),
      x = x,
      x_parents = x_parents,
      node_values = numeric(nrow(nodes)),
      node_values_map = seq(nrow(nodes)),
      nodes = nodes,
      node_graph = node_graph,
      variance = variance,
      parameters = parameters,
      parameter_map = seq_along(parameters),
      correlation_function = correlation_function,
      x_covariance = x_covariance,
      node_covariance = node_covariance
    )
  )
}

#' @describeIn update_spline Update a spline's covariance matrix
#' 
#' @param variance The new marginal variance of the spline. If missing, will
#'   use the pre-existing variance stored in the spline.
#' @param parameters The new parameters of the spline. If missing, will
#'   use the pre-existing parameters stored in the spline.
#' 
#' @export
update_spline_covariance<- function(
    spline,
    variance,
    parameters,
    ...
  ) {
  if( !missing(variance) ) spline$variance<- variance
  if( !missing(parameters) ) spline$parameters<- parameters
  mode<- "numeric"
  if( "advector" %in% class(spline$parameters) && requireNamespace("RTMB") ) mode<- "advector"
  if( mode == "advector" ) {
    spline$x_covariance<- as(spline$x_covariance, "adsparse")
    spline$node_covariance<- as(spline$node_covariance, "adsparse")
  } else {}

  x_m<- cbind(
    spline$x_covariance@i + 1,
    ip_to_j(spline$x_covariance@i, spline$x_covariance@p) + 1,
    0
  )
  x_m[, 3]<- spline$variance * apply(
    x_m,
    MARGIN = 1,
    function(row) spline$correlation_function(
      spline$x[row[[1]], ],
      spline$nodes[row[[2]], ],
      spline$parameters
    )
  )
  spline$x_covariance@x<- x_m[, 3]

  node_m<- cbind(
    spline$node_covariance@i + 1,
    ip_to_j(spline$node_covariance@i, spline$node_covariance@p) + 1,
    0
  )
  node_m[, 3]<- spline$variance * apply(
    node_m,
    MARGIN = 1,
    function(row) spline$correlation_function(
      spline$nodes[row[[1]], ],
      spline$nodes[row[[2]], ],
      spline$parameters
    )
  )
  spline$node_covariance@x<- node_m[, 3]

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
    spline$x_parents[node_x, ],
    MARGIN = 1,
    function(parents) which(parents == -1)
  )
  spline$values[node_x]<- spline$node_values[node_x_parents]

  # Update x values that are between nodes
  no_node_x<- which(!x_type)
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
      Sigma[1, 1]<- spline$variance
      Sigma[1, seq_along(p) + 1]<- spline$x_covariance[i, p]
      Sigma[seq_along(p) + 1, 1]<- 0
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
  if( "advector" %in% class(x) && requireNamespace("RTMB") ) mode<- "advector"
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
      cbind(spline$node_values[c(i, p)]),
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

#' Evaluate the penalized log-likelihood of a spline
#'
#' @param x The node values
#' @param spline The spline object
#' @param penalty How much should the variance be penalized?
#' @param log If TRUE, returns the log-likelihood.
#'
#' @return The log-likelihood of the spline evaluated at x. If x is a simref
#'   object, then the values of x will be replaced with simulated values.
#'
#' @export
dpspline<- function(
    x,
    spline,
    penalty = 0.01,
    log = TRUE
  ) {
  mode<- "numeric"
  if( "advector" %in% class(x) && requireNamespace("RTMB") ) mode<- "advector"
  ll<- 0
  if( mode == "advector" ) {
    ll<- ll + RTMB::dexp(
      spline$variance,
      penalty,
      TRUE
    )
  } else {
    ll<- ll + stats::dexp(
      spline$variance,
      penalty,
      TRUE
    )
  }
  ll<- ll + dspline(x, spline, log = TRUE)
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
#' 
#' @return The value(s) of the spline at x
#'
#' @export
nns<- function(
    x,
    spline
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
  return(spline$values[idx])
}

