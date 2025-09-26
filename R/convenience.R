#' Convert a distance matrix to a directed acyclic graph
#'
#' @param d 
#'     A distance matrix, or of class "dist"
#' @param k 
#'     The number of parents for each node
#'
#' @return
#'     A directed acyclic graph from the igraph package
distance_matrix_to_dag<- function(
        d,
        k = 4
    ) {
    if( "dist" %in% class(d) ) d<- as.matrix(d)
    order<- rep(NA, nrow(d))
    min_d<- numeric(nrow(d))
    order[1]<- 1
    min_d[1]<- Inf
    min_d[-1]<- d[1, -1]
  
    i<- 2
    while( any(is.na(order)) ) {
        order[i]<- which.min(min_d)
        min_d[order[i]]<- Inf
        new_d<- d[order[i], -order[seq(i)]]
        replace_idx<- (new_d < min_d[-order[seq(i)]])
        min_d[-order[seq(i)]][replace_idx]<- new_d[replace_idx]
        i<- i + 1
    }

    parent_list<- vector(mode = "list", nrow(d))
    parent_list[order]<- seq(nrow(d)) |>
        lapply(
            function(i) {
                if( i == 1 ) return(numeric(0))
                k<- min(k, i - 1)
                parent_d<- d[order[i], order[seq_len(i - 1)]]
                partial_sort<- parent_d |>sort(partial = k)
                parents<- (parent_d <= partial_sort[k]) |> which()
                return( order[parents] )
            }
        )

    m<- matrix(0, nrow = nrow(d), ncol = nrow(d))
    for( i in seq(nrow(d)) ) m[parent_list[[i]], i]<- 1
    dag<- m |> igraph::graph_from_adjacency_matrix(mode = "directed")
    return(dag)
}

#' Find each pair of nodes whose correlation needs to be computed.
#'
#' @param dag 
#'     An directed acyclic graph as an igraph object
#'
#' @return 
#'     A data.frame giving the list of every node that appears in the same
#'     parent list or are parent/children.
get_pairs<- function(
        dag
    ) {
    pairs<- seq_along(igraph::V(dag)) |>
        lapply(
            function(v) {
                parents<- dag |> igraph::neighbors(v, "in")
                pairs<- matrix(v, ncol = 2)
                if( length(parents) > 0 ) {
                    pairs<- rbind(
                        pairs,
                        c(v, parents) |> sort() |> combn(2) |> t()
                    )
                }
                colnames(pairs)<- c("i", "j")
                return(pairs)
            }
        ) |>
        do.call(rbind, args = _) |>
        unique() |>
        as.matrix()
    return(pairs)
}



#' Conditional mean and variance from joint Gaussian variables
#'
#' @param x 
#'     An n x k matrix. Columns index independent replicates.
#' @param joint_mean 
#'     An n x k mean matrix. Columns index independent replicates.
#' @param joint_covariance 
#'     An n x n covariance matrix.
#' @param predicted 
#'     Indices declaring which of the n elements should be predicted.
#' @param predictors 
#'     Indices declaring which of the n elements should be conditioned on.
#'
#' @return 
#'     A list with the conditional mean of variance of 
#'     x\[predicted, \] given x\[predictors, \]
#' 
#' @export
conditional_normal<- function(
        x,
        joint_mean,
        joint_covariance,
        predicted = 1,
        predictors = x |> nrow() |> seq() |> tail(-1)
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
                covariance = joint_covariance[
                    predicted, 
                    predicted, 
                    drop = FALSE
                ]
            )
        )
    }

    SigmaAA<- joint_covariance[predicted, predicted, drop = FALSE]
    SigmaBA<- joint_covariance[predictors, predicted, drop = FALSE]
    SigmaBB<- joint_covariance[predictors, predictors, drop = FALSE]
    SigmaAB_BBinv<- solve(SigmaBB, SigmaBA) |> t()

    return(
        list(
            mean = joint_mean[predicted, , drop = FALSE] + 
                SigmaAB_BBinv %*% (x - joint_mean)[predictors, , drop = FALSE],
            covariance = SigmaAA - SigmaAB_BBinv %*% SigmaBA
        )
    )
}
