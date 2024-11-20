adaptive_nodes<- function(
        data,
        k = 10
    ) {
    if( !("matrix" %in% class(data)) ) data<- as.matrix(data)
    bandwidth<- cov(data) * nrow(data) ^ (-1/(ncol(data) + 4))
    precision<- solve(bandwidth)
    kernel<- function(x, bandwidth) RTMB::dmvnorm(x, Sigma = bandwidth)
    kde<- function(x) {
        if( !("matrix" %in% class(x)) ) x<- as.matrix(x)
        ans<- apply(
            x,
            MARGIN = 1,
            function(y) {
                centered_nodes<- sweep(
                    data,
                    MARGIN = 2,
                    STATS = y
                )
                return( mean(kernel(centered_nodes, precision)) )
            }
        )
        return(ans)
    }
    kde_potential_energy<- function(nodes) -sum(kde(nodes))

    node_distance<- function(nodes, prec) {
        pairs<- t(combn(nrow(nodes), 2))
        ans<- apply(
            pairs,
            MARGIN = 1,
            function(ij) {
                diff<- nodes[ij[[1]], ] - nodes[ij[[2]], ]
                return(.Machine$double.eps + sqrt(t(diff) %*% precision %*% diff ))
            }
        )
        return(ans)
    }
    repulsion_potential_energy<- function(nodes) {
        const * sum(1 / node_distance(nodes))
    }
    potential_energy<- function(nodes) kde_potential_energy(nodes) + repulsion_potential_energy(nodes)

    opt<- nlminb(
        c(data[sample(nrow(data), size = k), ]),
        objective = function(x) {
            nodes<- matrix(x, nrow = k, ncol = ncol(data))
            return(potential_energy(nodes))
        }
    )
}


data<- MASS::mvrnorm(
    n = 100,
    mu = cbind(c(0, 0)),
    Sigma = matrix(c(2, 0.7, 0.7, 1), nrow = 2)
)