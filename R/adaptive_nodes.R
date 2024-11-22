library(profvis)

adaptive_nodes<- function(
        data,
        n_nodes = 10,
        max.it = 100
    ) {
    if( !("matrix" %in% class(data)) ) data<- as.matrix(data)
    cov<- cov(data)
    precision<- solve(cov)
    data_precision<- precision * nrow(data) ^ (1 / (ncol(data) + 4))
    node_precision<- precision * n_nodes ^ (1 / (ncol(data) + 4))
    data_cov<- solve(data_precision)
    node_cov<- solve(node_precision)
    dfun<- RTMB::MakeTape(
        function(args) {
            n<- sqrt(1 + length(args)) - 1
            x<- args[1:n]
            y<- args[1:n + n]
            precision<- RTMB::matrix(
                args[-(1:(2 * n))],
                nrow = n
            )
            d<- sqrt(
                t(x - y) %*% precision %*% (x - y)
            )
            return(d[1, 1])
        },
        c(data[1, ], data[2, ], node_precision)
    )
    distance_function<- function(x, y, precision) {
        d<- dfun(c(x, y, precision))
        return(.Machine$double.eps + d)
    }
    data_distance<- function(nodes) {
        pairs<- expand.grid(
            i = seq(nrow(nodes)),
            j = seq(nrow(data))
        )
        ans<- apply(
            pairs,
            MARGIN = 1,
            function(ij) distance_function(
                nodes[ij[[1]], ],
                data[ij[[2]], ],
                data_precision
            )
        )
        d<- matrix(0, nrow = nrow(nodes), ncol = nrow(data))
        d[as.matrix(pairs)]<- ans
        return(d)
    }
    data_penalty<- function(nodes) {
        # Want to penalize data points that are far from their nearest node and nodes that are far from the nearest data point
        d<- data_distance(nodes)
        min_node_d<- apply(d, MARGIN = 1, min)
        min_data_d<- apply(d, MARGIN = 2, min)
        return(sum(min_node_d) + sum(min_data_d) + 10 * max(min_node_d))
    }
    node_distance<- function(nodes) {
        if( nrow(nodes) == 1 ) return(numeric(0))
        pairs<- t(combn(nrow(nodes), 2))
        ans<- apply(
            pairs,
            MARGIN = 1,
            function(ij) distance_function(
                nodes[ij[[1]], ],
                nodes[ij[[2]], ],
                node_precision
            )
        )
        m<- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))
        m[pairs]<- ans
        return(m)
    }
    node_penalty<- function(nodes) {
        # Want to penalize nodes that are close together
        if( nrow(nodes) <= 1 ) return(0)
        d<- node_distance(nodes)
        d[lower.tri(d, diag = TRUE)]<- NA
        min_d<- apply(d[, -1, drop = FALSE], MARGIN = 2, min, na.rm = TRUE)
        return(sum(1 / min_d))
    }
    total_penalty<- function(nodes) data_penalty(nodes) + node_penalty(nodes)

    penalty_history<- numeric(max.it + 1)
    nodes<- data[sample(nrow(data), size = n_nodes), , drop = FALSE]
    penalty<- total_penalty(nodes)
    penalty_history[1]<- penalty
    for( i in seq(max.it) ) {
        temperature<- 0.5 * ((max.it - i) / max.it)^4
        type<- sample(
            c("Jitter", "Replace"),
            size = sample(2, size = 1, prob = c(0.8, 0.2)),
            prob = c(
                0.8,
                0.2
            )
        )
        if( "Jitter" %in% type ) {
            new_nodes<- nodes + MASS::mvrnorm(
                n = nrow(nodes),
                mu = cbind(c(0, 0)),
                Sigma = 0.001 * node_cov
            )
        }
        if( "Replace" %in% type ) {
            new_nodes<- nodes
            new_nodes[
                sample(nrow(new_nodes), size = 1),

            ]<- data[
                sample(nrow(data), size = 1),

            ]
        }
        if( "Add" %in% type ) {
            new_nodes<- rbind(
                nodes,
                data[sample(nrow(data), size = 1), , drop = FALSE]
            )
        }
        new_penalty<- total_penalty(new_nodes)
        acceptance<- exp(-(new_penalty - penalty) / temperature)
        if( (new_penalty < penalty) | (runif(1) < acceptance) ) {
            nodes<- new_nodes
            penalty<- new_penalty
        }
        penalty_history[i + 1]<- penalty
    }
    return(list(nodes = nodes, penalty_history = penalty_history))
}


data<- rbind(
    MASS::mvrnorm(
        n = 100,
        mu = cbind(c(4, 4)),
        Sigma = matrix(c(2, -1.7, -1.7, 2), nrow = 2)
    ),
    MASS::mvrnorm(
        n = 50,
        mu = cbind(c(0, 0)),
        Sigma = matrix(c(2, 0.7, 0.7, 1), nrow = 2)
    ),
    MASS::mvrnorm(
        n = 50,
        mu = cbind(c(2, 2)),
        Sigma = matrix(c(0.2, 0, 0, 0.2), nrow = 2)
    ),
    MASS::mvrnorm(
        n = 50,
        mu = cbind(c(-3, 5)),
        Sigma = matrix(c(3, 1, 1, 0.4), nrow = 2)
    ),
    MASS::mvrnorm(
        n = 1,
        mu = cbind(c(10, 10)),
        Sigma = matrix(c(1, 0, 0, 1), nrow = 2)
    )
)

profvis({
    nodes<- adaptive_nodes(data, n_nodes = 30, max.it = 1000)
})


plot(data, ylim = c(-6, 12), xlim = c(-6, 12), pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8)
points(nodes$nodes, pch = 19, col = 2, cex = 0.8)






