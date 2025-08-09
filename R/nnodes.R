#' Find a good set of nodes to cover a dataset
#' 
#' @param data 
#'     A matrix of data points. Each row should be a data point.
#' @param n_nodes 
#'     The number of nodes to cover the dataset.
#' @param max.it 
#'     The number of iterations to use when finding the best set of nodes.
#' @param jitter_magnitude 
#'     The relative step size for randomly sampled jitters
#' @param create_graph 
#'     Should a directed acyclic graph between the nodes be returned?
#' @param n_parents 
#'     If create_graph is TRUE, the number of parents for each node.
#' 
#' @return 
#'     A list giving the nodes, an LT matrix for computing mahalanobis distance,
#'     the objective function trajectory, and possibly the node graph.
#' 
#' @export
nnodes<- function(
        data,
        n_nodes = 10,
        max.it = 100,
        jitter_magnitude = 0.0001,
        create_graph = FALSE,
        n_parents = 4
    ) {
    if( !requireNamespace("RANN", quietly = TRUE) ) {
        stop("Must have the RANN package installed to use the annodes function.")
    }
    if( !("matrix" %in% class(data)) ) data<- data |> as.matrix()
    data<- data[data |> complete.cases(), , drop = FALSE]

    cov<- data |> cov()
    LT<- cov |> symsqrt(invert = TRUE)
    trans_data<- data %*% LT
    find_nearest_distances<- function(nodes) {
        trans_nodes<- nodes %*% LT
        data_to_node_d<- RANN::nn2(
            trans_nodes,
            trans_data,
            k = 1
        )$nn.dists[, 1]
        node_to_data_d<- RANN::nn2(
            trans_data,
            trans_nodes,
            k = 1
        )$nn.dists[, 1]
        node_to_node_d<- RANN::nn2(
            trans_nodes,
            trans_nodes,
            k = 2
        )$nn.dists[, -1]
        return(
            list(
                data_to_node = data_to_node_d,
                node_to_data = node_to_data_d,
                node_to_node = node_to_node_d
            )
        )
    }
    penalty_function<- function(dlist) {
        mean(dlist$data_to_node) + 
            mean(dlist$node_to_data) - 
            mean(dlist$node_to_node)
    }

    penalty_history<- data.frame(
        penalty = numeric(max.it + 2),
        operation = ""
    )
    nodes<- data[data |> nrow() |> sample(size = n_nodes), , drop = FALSE]
    dlist<- nodes |> find_nearest_distances()
    penalty<- dlist |> penalty_function()
    best_nodes<- nodes
    best_penalty<- penalty
    penalty_history$penalty[1]<- penalty
    penalty_history$operation[1]<- "Start"
    for( i in seq(max.it) ) {
        T<- temperature(i, max.it, 1, "decay")
        operation<- sample(
            c("Jitter", "Replace"),
            size = sample(2, size = 1, prob = c(0.8, 0.2)),
            prob = c(
                0.8,
                0.2
            )
        )
        new_nodes<- nodes
        if( "Jitter" %in% operation ) {
            new_nodes<- nodes |>
                jitter(
                    jitter_magnitude,
                    cov
                )
        }
        if( "Replace" %in% operation ) {
            new_nodes<- nodes |>
                replace(
                    data,
                    dlist
                )
        }
        new_dlist<- new_nodes |> find_nearest_distances()
        new_penalty<- new_dlist |> penalty_function()
        if( accept(new_penalty, penalty, T, 0.1) ) {
            nodes<- new_nodes
            penalty<- new_penalty
        } else {
            operation<- "Keep"
        }
        if( penalty < best_penalty ) {
            best_nodes<- nodes
            best_penalty<- penalty
        } else if( restart(new_penalty, best_penalty, T, 0.3) ) {
            nodes<- best_nodes
            penalty<- best_penalty
            operation<- "Restart"
        }
        penalty_history$penalty[i + 1]<- penalty
        penalty_history$operation[i + 1]<- operation |> paste(collapse = "")
    }
    nodes<- best_nodes
    penalty_history$penalty[i + 2]<- best_penalty
    penalty_history$operation[i + 2]<- "Best"

    ans<- list(
        nodes = best_nodes,
        LT = LT,
        penalty = penalty_history
    )
    if( create_graph ) {
        ans$graph<- (nodes %*% LT) |> 
            dist() |> 
            distance_matrix_to_dag(
                k = n_parents
            )
    }
    return( ans )
}

#' Compute symmetric square roots of positive definite matrices
#' 
#' @param m 
#'     A positive definite matrix
#' @param invert 
#'     If TRUE, the take the symmetric square root of m^-1
#' 
#' @return 
#'     A symmetric matrix x where x %*% x == m (or m^-1)
#' 
#' @export
symsqrt<- function(
        m,
        invert = FALSE
    ) {
    eig<- m |> eigen()
    eig$values<- eig$values^(0.5 * ifelse(invert, -1, 1))
    val_vecinv<- solve(
            a = eig$vector |> t(),
            b = diag(eig$values, nrow = nrow(m)) |> t()
        ) |>
        t()
    ans<- eig$vectors %*% val_vecinv
    return(ans)
}

temperature<- function(
        i,
        max.it,
        T0 = 1,
        type = c(
            "default",
            "decay",
            "exp",
            "fast",
            "boltzmann"
        )
    ) {
    T<- switch(
        type[[1]],
        "default" = T0 * 0.5 * (1 - (i - 1) / max.it)^4,
        "decay" = T0 * (1 - (i - 1) / max.it),
        "exp" = T0 * 0.95^i,
        "fast" = T0 * (1 / i),
        "boltzmann" = T0 * (1 / log(i))
    )
    return(T)
}

accept<- function(
        new_penalty,
        old_penalty,
        T,
        max_probability
    ) {
    if( new_penalty < old_penalty ) return(TRUE)
    p<- exp( -(new_penalty - old_penalty) * (1 / T))
    p<- max_probability * p
    return(sample(c(TRUE, FALSE), 1, prob = c(p, 1 - p)))
}

restart<- function(
        new_penalty,
        best_penalty,
        T,
        max_probability
    ) {
    p<- 1 - exp( -(new_penalty - best_penalty) * T)
    p<- max_probability * p
    return(sample(c(TRUE, FALSE), 1, prob = c(p, 1 - p)))
}

jitter<- function(
        nodes,
        magnitude,
        covariance
    ) {
    nodes<- nodes + MASS::mvrnorm(
        n = nrow(nodes),
        mu = cbind(numeric(ncol(nodes))),
        Sigma = magnitude * covariance
    )
    return(nodes)
}

replace<- function(
        nodes,
        data,
        dlist
    ) {
    node_penalty<- dlist$node_to_data - dlist$node_to_node
    node_penalty<- node_penalty + min(node_penalty)
    data_penalty<- dlist$data_to_node
    delete_idx<- sample(
        nrow(nodes),
        size = 1,
        prob = node_penalty / sum(node_penalty)
    )
    add_idx<- sample(
        nrow(data),
        size = 1,
        prob = data_penalty / sum(data_penalty)
    )
    nodes[delete_idx, ]<- data[add_idx, ]
    return(nodes)
}