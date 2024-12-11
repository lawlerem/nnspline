#' Convert an (i, p) sparse matrix index to an (i, j) index
#' 
#' @param i the i index
#' @param p the p index
#' 
#' @return A vector giving the column index of the entries
ip_to_j<- function(
    i,
    p
  ) {
  j<- findInterval(
    seq(i) - 1,
    p[-1]
  )
  return(j)
}

#' Convert a sparse matrix to a dense matrix, preserving RTMB evaluation mode
#'
#' @param m A sparse matrix
#' @param i An integer vector to index rows to subset
#' @param j An integer vector to index columns to subset
#' @param m_ij Optional. An integer vector giving the vector position of the non-zero values of m if treated as a dense matrix.
#'
#' @return A dense matrix with the appropriate evaluation class
adsparse_to_admatrix<- function(m, i, j, m_ij) {
  if( !requireNamespace("RTMB") ) stop("must have RTMB package to use this function.")
  if( missing(i) ) i<- seq_len(m@Dim[[1]])
  if( missing(j) ) j<- seq_len(m@Dim[[2]])
  if( any(i > m@Dim[[1]]) ) stop("i is greater than the number of rows.")
  if( any(j > m@Dim[[2]]) ) stop("j is greater than the number of columns.")

  ans<- RTMB::AD(numeric(length(i) * length(j)))
  dim(ans)<- c(length(i), length(j))
  if( missing(m_ij) ) m_ij<- m@i + ip_to_j(m@i, m@p) * m@Dim[[1]]
  mx<- m@x
  
  ans_ij<- c(outer(i - 1, j - 1, function(X, Y) X + Y * m@Dim[1]))
  ans_ind<- cbind(
    rep_len(seq_along(i), length(i) * length(j)),
    rep(seq_along(j), each = length(i))
  )
  matches<- match(ans_ij, m_ij)
  isgood_matches<- !is.na(matches)
  ans[ans_ind[isgood_matches, ]]<- mx[matches[isgood_matches]]

  return(ans)
}


#' Find the joint covariance matrix for an x/node and its parents
#' 
#' @param spline The spline
#' @param idx The index of either an x or node
#' @param parent_idx The index of the parents of idx
#' @param prediction_type Either "x" if idx refers to an x or "node" if idx refers to a node
#' @param mode Either "numeric" or "advector"
#' 
#' @return A dense covariance matrix
get_joint_covariance<- function(
    spline, 
    idx,
    parent_idx,
    prediction_type = "node",
    mode = "numeric"
  ) {
  if( prediction_type == "node" ) {
    if( mode == "numeric" ) {
      Sigma<- as.matrix(spline$node_covariance[c(idx, parent_idx), c(idx, parent_idx)])
    } else if( mode == "advector" ) {
      Sigma<- adsparse_to_admatrix(
        m = spline$node_covariance,
        i = c(idx, parent_idx),
        j = c(idx, parent_idx),
        m_ij = spline$node_ij
      )
    }
  } else if( prediction_type == "x" ) {
    if( mode == "numeric" ) {
      Sigma<- as.matrix(spline$node_covariance[c(1, parent_idx), c(1, parent_idx)])
      Sigma[1, ]<- as.matrix(spline$x_covariance[idx, c(1, parent_idx + 1)])
    } else if( mode == "advector" ) {
      Sigma<- adsparse_to_admatrix(
        spline$node_covariance,
        c(1, parent_idx),
        c(1, parent_idx),
        spline$node_ij
      )
      Sigma[1, ]<- adsparse_to_admatrix(
        m = spline$x_covariance,
        i = idx,
        j = c(1, parent_idx + 1),
        m_ij = spline$x_ij
      )
    }
  }
  Sigma[-1, 1]<- 0
  Sigma<- Sigma + t(Sigma)
  diag_ind<- cbind(seq(nrow(Sigma)), seq(nrow(Sigma)))
  Sigma[diag_ind]<- Sigma[diag_ind] / 2
  
  return(Sigma)
}

#' Convert a distance matrix to a directed acyclic graph
#'
#' @param d A distance matrix, or of class "dist"new_
#' @param k The number of parents for each node
#'
#' @return A directed acyclic graph from the igraph package
distance_matrix_to_dag<- function(d, k = 4) {
  if( "dist" %in% class(d) ) d<- as.matrix(d)
  dag<- .distance_matrix_to_dag(d, k) # from Rcpp nnspline.cpp
  order<- dag$order
  parents<- dag$parent_list
  m<- matrix(
    0,
    nrow = length(order),
    ncol = length(order)
  )
  for( i in seq_along(parents) ) {
    m[parents[[i]], i]<- 1
  }
  inverse_order<- seq_along(order)[order(order)]
  m<- m[inverse_order, inverse_order]
  dag<- igraph::graph_from_adjacency_matrix(m, mode = "directed")
  return(dag)
}

#' Find each pair of nodes whose correlation needs to be computed.
#'
#' @param dag An directed acyclic graph as an igraph object
#' @param x_parents A sparse matrix which has a 1 in entry (i, j) if
#'   node i is a parent of x\[j\]
#'
#' @return A data.frame giving the list of every node
#'   that appear in the same parent list or are parent/children.
get_pairs<- function(dag, x_parents) {
  dag_pairs<- lapply(
    seq_along(igraph::V(dag)),
    function(v) {
      parents<- igraph::neighbors(dag, v, "in")
      pairs<- matrix(v, ncol = 2)
      if( length(parents) > 0 ) {
        pairs<- rbind(
          pairs,
          t(combn(sort(c(v, parents)), 2))
        )
      } else {}
      pairs<- as.data.frame(pairs)
      colnames(pairs)<- c("i", "j")
      return(pairs)
    }
  )
  x_pairs<- apply(
    x_parents,
    MARGIN = 1,
    function(row) {
      parents<- which(row > 0)
      if( length(parents) == 0 ) {
        return(
          data.frame(
            i = integer(0),
            j = integer(0)
          )
        )
      } else {}
      pairs<- t(combn(sort(parents), 2))
      pairs<- as.data.frame(pairs)
      colnames(pairs)<- c("i", "j")
      return(pairs)
    },
    simplify = FALSE
  )
  pairs<- do.call(
    rbind,
    c(
      dag_pairs,
      x_pairs
    )
  )
  pairs<- unique(pairs)
  return(pairs)
}



#' Conditional mean and variance from joint Gaussian variables
#'
#' @param x An n x k matrix. Columns index independent replicates.
#' @param joint_mean An n x k mean matrix. Columns index independent replicates.
#' @param joint_covariance An n x n covariance matrix. Can be a sparse matrix.
#' @param predicted Indices declaring which of the n elements should be predicted.
#' @param predictors Indices declaring which of the n elements should be conditioned on.
#'
#' @return A list with the conditional mean of variance of x\[predicted, \] given x\[predictors, \]
#' 
#' @export
conditional_normal<- function(
    x,
    joint_mean,
    joint_covariance,
    predicted = 1,
    predictors = tail(seq(nrow(x)), -1)
) {
  mode<- "numeric"
  if( "advector" %in% class(joint_covariance) && requireNamespace("RTMB") ) mode<- "advector"
  if( length(predictors) == 0 ) {
    return(
      list(
        mean = joint_mean[predicted, , drop = FALSE],
        covariance = joint_covariance[predicted, predicted, drop = FALSE]
      )
    )
  } else {}

  SigmaAA<- joint_covariance[predicted, predicted, drop = FALSE]
  SigmaAB<- t(joint_covariance[predictors, predicted, drop = FALSE])
  SigmaBB<- joint_covariance[predictors, predictors, drop = FALSE]

  if( mode == "advector" ) {
    SigmaAB_BBinv<- SigmaAB %*% RTMB::solve(SigmaBB)
  } else {
    SigmaAB_BBinv<- SigmaAB %*% solve(SigmaBB)
  }
  return(
    list(
      mean = joint_mean[predicted, , drop = FALSE] + SigmaAB_BBinv %*% (x[predictors, , drop = FALSE] - joint_mean[predictors, , drop = FALSE]),
      covariance = SigmaAA - SigmaAB_BBinv %*% t(SigmaAB)
    )
  )
}




