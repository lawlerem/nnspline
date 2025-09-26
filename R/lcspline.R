#' Create a lcspline object
#' 
#' @param x 
#'     A numeric matrix giving the locations (as rows) for which the spline
#'     values are desired
#' @param n_parents 
#'     The number of parents each node should have.
#' @param correlation_function 
#'     A function used to compute correlations between spline inputs. 
#'     Should have the arguments x1 and x2, vectors giving the spline inputs, 
#'     and a stretch parameter equivalent to the traditional range parameter of
#'     Matern covariance functions.
#' @param stretch_limits
#'     Optional. Lower and upper limits of stretch_sequence. If missing, a
#'     lower limit of 0 is used and a sensible upper limit is found based on
#'     a random subset of the spline neighbourhood structure.
#' @param n_stretch
#'     The number of stretch values used. Values will be equally spaced between
#'      min(stretch_limits) and max(stretch_limits).
#' @param stretch_sequence
#'     Optional. A numeric vector giving a sequence of stretch values used to 
#'     precompute the spline structure. This value overrides both stretch_limits
#'     and n_stretch.
#' @param graph 
#'     Optional. If supplied, it needs to be a directed acyclic graph from the
#'     igraph package. The graph described the structure of the spline by
#'     defining a sequence of conditional factorizations of the joint likelihood.
#' @param ...
#'     Arguments to pass to find_upper_stretch.
#' 
#' @return
#'     A lcspline object as a list with the following elements:
#'   * map 
#'       Mapping indices for the node values. See ?TMB::MakeADFun.
#'   * x 
#'       The locations at which the spline is evaluated.
#'   * graph 
#'       An igraph directed acyclic graph represnting the node parent structure.
#'   * n_parents 
#'       The number of parents for each location
#'   * cond_mean
#'       A list of sparse matrices giving the conditional mean structure of the
#'       spline for each stretch value.
#'   * cond_var
#'       A list giving the conditional variances of each node for each stretch
#'       value.
#'   * stretch_sequence
#'       The sequence of range parameters used to compute LC.
#'   * mixtape
#'       An AD-compatible function used to compute mixture probabilities used
#'       to interpolate between precomputed values of LC.
#'   * correlation_function
#'       The function used to compute correlations between input points.
#' 
#' @export
create_lcspline<- function(
        x,
        n_parents = 4,
        correlation_function = function(x1, x2, stretch) {
            d<- sqrt(sum((x1 - x2)^2))
            if( stretch == 0 ) return( (d == 0) * 1 )
            (1 + d/stretch) * exp( -d/stretch )
        },
        stretch_limits,
        n_stretch = 5,
        stretch_sequence,
        graph,
        ...
    ) {
    if( !("matrix" %in% class(x)) ) x<- x |> matrix(ncol = 1)
    n_parents<- min(n_parents, nrow(x))
    
    if( missing(graph) ) {
        graph<- x |> dist() |> distance_matrix_to_dag(k = n_parents)
    }
    if( missing(stretch_sequence) & missing(stretch_limits) ) {
        if( missing(stretch_limits) ) {
            upper_stretch<- find_upper_stretch(
                x,
                graph,
                correlation_function,
                ...
            )
            stretch_limits<- c(0, upper_stretch)
        }
        stretch_sequence<- seq(
            min(stretch_limits), 
            max(stretch_limits), 
            length.out = n_stretch
        )
    }
    stretch_sequence<- stretch_sequence |> sort()

    spline<- list(
        map = x |> nrow() |> seq(),
        x = x,
        graph = graph,
        n_parents = n_parents,
        cond_mean = list(),
        cond_var = list(),
        stretch_sequence = 0,
        mixtape = function() {},
        correlation_function = correlation_function
    )
    spline<- spline |> cache_lc(stretch_sequence)
    spline
}

mixtape_recorder<- function(nodes) {
    if( requireNamespace("RTMB", quitely = TRUE) ) {
        mixtape<- RTMB::MakeTape(
            \(t) {
                nodes<- nodes |> RTMB::AD()
                t<- (max(nodes) - min(nodes)) * t + min(nodes)
                p<- ((tail(nodes, -1) - t) / (tail(nodes, -1) - head(nodes, -1))) |>
                    RTMB::sapply(\(x) min(c(x, 1))) |> 
                    RTMB::sapply(\(x) max(c(x, 0))) |>
                    c(1)
                pp<- p
                pp[1]<- p[1]
                for( i in pp |> seq_along() |> tail(-1) ) {
                    pp[i]<- p[i] - sum(pp[seq(i - 1)])
                }
                pp
            },
            c(0.5)
        )
    } else {
        mixtape<- function(t) {
            t<- (max(nodes) - min(nodes)) * t + min(nodes)
            p<- ((tail(nodes, -1) - t) / (tail(nodes, -1) - head(nodes, -1))) |>
                sapply(\(x) min(c(x, 1))) |> 
                sapply(\(x) max(c(x, 0))) |>
                c(1)
            pp<- p
            pp[1]<- p[1]
            for( i in pp |> seq_along() |> tail(-1) ) {
                pp[i]<- p[i] - sum(pp[seq(i - 1)])
            }
            pp
        }
    }
    return(mixtape)
}

#' Find a sensible upper limit for the range parameter
#' 
#' When assuming each variable has a marginal variance of 1, the conditional
#'     variance of each variable given its parents is bounded between 0 and 1,
#'     and decreases from 1 to 0 as the range parameter increases as long as
#'     the correlation function is monotonic decreasing.
#' When the conditional variance is sufficiently small, increasing the range
#'     parameter has a negligble effect on the conditional distribution and
#'     can lead to nearly flat neighbourhoods of the likelihood function.
#' This function find a sensible upper limit for the range parameter by finding
#'     a value p such that the conditional variance for the variable is 
#'     sufficiently small for an accurate conditional distribution but not small
#'     enough for a flat likelihood.
#'
#' @param x The spline x coordinates
#' @param graph The spline graph
#' @param correlation_function The correlation function
#' @param ignore_threshold If x has this many parents or fewer, it will not be
#'     used when finding the upper limit.
#' @param subsample An integer giving the number of randomly sampled nodes used
#'     to find the upper limit.
#' @param target_variance The target value for the conditional variance at the
#'     the upper limit of the range parameter.
#' @param quantile The quantile used to pick the aggregate upper range from the
#'     the individual range estimates.
#' 
#' @return A numeric value giving a sensible upper limit for the range parameter.
#' 
#' @export
find_upper_stretch<- function(
        x,
        graph,
        correlation_function,
        ignore_threshold = igraph::max_degree(graph, mode = "in") - 1,
        subsample = 30,
        target_variance = 0.01,
        quantile = 0.8
    ) {
        ignore_threshold<- min(
            ignore_threshold,
            graph |> igraph::max_degree(mode = "in") |> (\(x) x - 1)()
        )
    is_ignored<- graph |> 
        igraph::degree(mode = "in") |> 
        (\(degree) degree <= ignore_threshold)()
    target_x<- x |> 
        nrow() |> 
        seq() |>
        _[!is_ignored] |>
        (\(i) sample(i, min(subsample, length(i))))()
    obj<- target_x |>
        lapply(
            \(v) {
                p<- graph |> igraph::neighbors(v, "in")
                v<- c(v, p)
                n<- length(v)
                min_d<- x[p, ] |> dist() |> as.matrix() |> _[-1, 1] |> min()

                fn<- function(range) {
                    S<- matrix(0, n, n)
                    for( i in seq(n) ) {
                        for( j in seq(i) ) {
                            S[i, j]<- S[j, i]<- correlation_function(
                                x[v[i], ],
                                x[v[j], ],
                                range
                            )
                        }
                    }
                    AA<- S[1, 1, drop = FALSE]
                    BA<- S[-1, 1, drop = FALSE]
                    BB<- S[-1, -1, drop = FALSE]
                    AB_BBinv<- t(solve(BB, BA))
                    var<- (AA - AB_BBinv %*% BA) |> as.numeric()
                    return( (var - target_variance)^2 )
                }
                return(
                    list(
                        par = min_d,
                        fn = fn
                    )
                )
            }
        )
    stretch<- obj |> 
        lapply(\(x) nlminb(x$par, x$fn)) |>
        lapply(`[[`, "par") |>
        do.call(c, args = _) |>
        quantile(probs = 0.8) |>
        unname()
    return(stretch)
}

#' Precompute spline structure
#' 
#' @param stretch_sequence A numeric vector giving the range parameters used to
#'     precompute the spline structure.
#' @param spline A lcspline object.
#' 
#' @return A copy of spline with an updated precomputed spline structure
#' 
#' @export
cache_lc<- function(spline, stretch_sequence) {
    spline$stretch_sequence<- stretch_sequence
    LC<- stretch_sequence |>
        lapply(
            compute_single_lc,
            spline$correlation_function,
            spline$graph,
            spline$x
        )
    spline$cond_mean<- LC |> lapply(`[[`, "plc")
    spline$cond_var<- LC |> lapply(`[[`, "var")
    spline$mixtape<- mixtape_recorder(stretch_sequence)
    return(spline)
}

compute_single_lc<- function(lc, corfun, graph, x) {
    p<- graph |> 
        igraph::V() |> 
        lapply(\(i) graph |> igraph::neighbors(i, "in") |> as.numeric())
    LC<- graph |>
        igraph::V() |>
        seq_along() |>
        lapply(
            \(idx) {
                v<- idx |> c(p[[idx]])
                n<- v |> length()
                S<- matrix(0, n, n)
                for( i in seq(n) ) {
                    for( j in seq(i) ) {
                        S[i, j]<- 
                            S[j, i]<- 
                            corfun(
                                x[v[i], ], 
                                x[v[j], ], 
                                lc
                            )
                    }
                }
                AA<- S[1, 1, drop = FALSE]
                BA<- S[-1, 1, drop = FALSE]
                BB<- S[-1, -1, drop = FALSE]
                if( n == 1 ) {
                    AB_BBinv<- matrix(0, nrow = 1, ncol = 0)
                } else {
                    AB_BBinv<- solve(BB, BA) |> t()
                }
                plc<- AB_BBinv |> 
                    as.numeric() |>
                    data.frame(
                        i = idx |> rep(length(p[[idx]])),
                        j = p[[idx]],
                        x = _
                    )
                var<- (AA - AB_BBinv %*% BA) |> as.numeric()
                list(
                    plc = plc, 
                    var = var
                )
            }
        )
    plc<- LC |>
        lapply(`[[`, "plc") |>
        do.call(rbind, args = _) |>
        subset(x != 0) |>
        (\(x) Matrix::sparseMatrix(
            i = x$i,
            j = x$j,
            x = x$x,
            dims = graph |> igraph::V() |> length() |> rep(2)
        ))()
    var<- LC |>
        lapply(`[[`, "var") |>
        do.call(c, args = _)
    return(list(plc = plc, var = var))
}

#' Evaluate the log-likelihood of a lcspline
#' 
#' @param x
#'     The output values of the spline.
#' @param spline
#'     A lcspline object.
#' @param stretch
#'     A numeric value between 0 and 1 determining the stretch of the spline.
#'     The value between 0 and 1 is shifted and scaled to the minimum and
#'     maximum values of the spline's stretch_sequence.
#' @param height
#'     A positive numeric value proportional to the height of the spline. The
#'     conditional variance is multiplied by height^2.
#' @param log
#'     If true, returns the log-likelihood.
#' 
#' @return
#'     The log-likelihood of the spline evaluated at x.
#' 
#' @export
dlcspline<- function(
        x, 
        spline, 
        stretch, 
        height, 
        log = TRUE
    ) {
    if( requireNamespace("RTMB", quietly = TRUE ) ) {
        dnorm<- RTMB::dnorm
    } else {
        dnorm<- stats::dnorm
    }

    ll<- 0
    mix<- spline$mixtape(stretch)
    cond_mean<- spline |>
        _$cond_mean |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)

    cond_var<- spline |>
        _$cond_var |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)
    if( x |> inherits("simref") ) {
        for( i in spline$graph |> igraph::topo_sort() |> as.numeric() ) {
            if( is.null(dim(x)) ) dim(x)<- c(length(x), 1)
            for( j in x |> ncol() |> seq() ) {
                p<- spline$graph |> igraph::neighbors(i, "in") |> as.numeric()
                ll<- ll + sum(
                    dnorm(
                        x[i, j],
                        sum(cond_mean[i, p] * x[p, j]),
                        height * sqrt(cond_var[i]),
                        log = TRUE
                    )
                )
            }
        }
    } else {
        mean<- (cond_mean %*% x) |> as.matrix()
        ll<- ll + sum(dnorm(x, mean, height * sqrt(cond_var), log = TRUE))
    }
    if( !log ) ll<- exp(ll)
    return( ll )
}

#' Simulated a lcspline
#' 
#' @param n
#'     The number of simulations to return.
#' @param spline
#'     The lcspline to simulate from.
#' @param stretch
#'     A numeric value between 0 and 1 determining the stretch of the spline.
#'     The value between 0 and 1 is shifted and scaled to the minimum and
#'     maximum values of the spline's stretch_sequence.
#' @param height
#'     A positive numeric value proportional to the height of the spline. The
#'     conditional variance is multiplied by height^2.
#' 
#' @return
#'     Values simulated from the spline.
#' 
#' @export
rlcspline<- function(n = 1, spline, stretch, height) {
    if( requireNamespace("RTMB", quietly = TRUE ) ) {
        dnorm<- RTMB::dnorm
    } else {
        dnorm<- stats::dnorm
    }
    values<- matrix(
        0,
        nrow = spline |> _$graph |> igraph::V() |> length(),
        ncol = n
    )
    mix<- spline$mixtape(stretch)
    cond_mean<- spline |>
        _$cond_mean |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)

    cond_var<- spline |>
        _$cond_var |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)
    for( i in spline$graph |> igraph::topo_sort() |> as.numeric() ) {
        for( j in n |> seq() ) {
            p<- spline$graph |> igraph::neighbors(i, "in") |> as.numeric()
            values[i, j]<- rnorm(
                    1,
                    sum(cond_mean[i, p] * values[p, j]),
                    height * sqrt(cond_var[i])
                )
        }
    }
    if( n == 1 ) values <- c(values)

    return( values )
}


