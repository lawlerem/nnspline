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
#' @param x_parents A list of sets of parent nodes
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
  x_pairs<- lapply(
    x_parents,
    function(p) {
      if( (length(p) == 1) && (p < 0) ) {
        return(
          data.frame(
            i = integer(0),
            j = integer(0)
          )
        )
      } else {}
      pairs<- t(combn(p, 2))
      pairs<- as.data.frame(pairs)
      colnames(pairs)<- c("i", "j")
      return(pairs)
    }
  )
  pairs<- do.call(
    rbind,
    c(
      dag_pairs,
      x_pairs
    )
  )
  pairs<- as.matrix(unique(pairs))
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
  if( "advector" %in% class(joint_covariance) && requireNamespace("RTMB") ) {
    mode<- "advector"
    solve<- RTMB::solve
  }
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
  SigmaAB_BBinv<- t(solve(SigmaBB, t(SigmaAB)))
  return(
    list(
      mean = joint_mean[predicted, , drop = FALSE] + SigmaAB_BBinv %*% (x[predictors, , drop = FALSE] - joint_mean[predictors, , drop = FALSE]),
      covariance = SigmaAA - SigmaAB_BBinv %*% t(SigmaAB)
    )
  )
}




