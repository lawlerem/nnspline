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
#'
#' @return A dense matrix with the appropriate evaluation class
adsparse_to_admatrix<- function(m) {
  if( !requireNamespace("RTMB") ) stop("must have RTMB package to use this function.")
  mm<- RTMB::AD(numeric(prod(m@Dim)))
  dim(mm)<- m@Dim
  i<- m@i
  j<- nnspline:::ip_to_j(m@i, m@p)
  ind<- cbind(i + 1, j + 1)
  x<- m@x
  mm[ind]<- x

  return(mm)
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
#' @param mu An n x k mean matrix. Columns index independent replicates.
#' @param Sigma An n x n covariance matrix. Can be a sparse matrix.
#' @param predicted Indices declaring which of the n elements should be predicted.
#' @param predictors Indices declaring which of the n elements should be conditioned on.
#'
#' @return A list with the conditional mean of variance of x\[predicted, \] given x\[predictors, \]
#' 
#' @export
conditional_normal<- function(
    x,
    mu,
    Sigma,
    predicted = 1,
    predictors = tail(seq(nrow(x)), -1)
) {
  mode<- "numeric"
  if( "advector" %in% class(Sigma) && requireNamespace("RTMB") ) mode<- "advector"
  if( length(predictors) == 0 ) {
    return(
      list(
        mean = mu[predicted, , drop = FALSE],
        covariance = Sigma[predicted, predicted, drop = FALSE]
      )
    )
  } else {}

  SigmaAA<- Sigma[predicted, predicted, drop = FALSE]
  SigmaAB<- t(Sigma[predictors, predicted, drop = FALSE])
  SigmaBB<- Sigma[predictors, predictors, drop = FALSE]

  if( mode == "advector" ) {
    SigmaAB_BBinv<- SigmaAB %*% RTMB::solve(SigmaBB)
  } else {
    SigmaAB_BBinv<- SigmaAB %*% solve(SigmaBB)
  }
  return(
    list(
      mean = mu[predicted, , drop = FALSE] + SigmaAB_BBinv %*% (x[predictors, , drop = FALSE] - mu[predictors, , drop = FALSE]),
      covariance = SigmaAA - SigmaAB_BBinv %*% t(SigmaAB)
    )
  )
}




