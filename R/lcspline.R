#' Create a lcspline object
#' 
#' @param x 
#'     A numeric matrix giving the locations (as rows) for which the spline
#'     values are desired
#' @param n_parents 
#'     The number of parents each node should have.
#' @param covariance_function 
#'     A function used to compute covariances between spline inputs. 
#'     Should have the arguments x1 and x2, vectors giving the spline inputs, 
#'     and a range parameter.
#' @param lc_sequence
#'     A numeric vector giving a sequence of range parameters used to precompute
#'     the spline structure.
#' @param graph 
#'     An optional argument. If supplied, it needs to be an igraph object 
#'     representing a directed acyclic graph.
#' 
#' @return
#'     A lcspline object as a list with the following elements:
#'   * values 
#'       The output values of the spline at x.
#'   * map 
#'       Mapping indices for the node values. See ?TMB::MakeADFun.
#'   * x 
#'       The locations at which the spline is evaluated.
#'   * graph 
#'       An igraph directed acyclic graph represnting the node parent structure.
#'   * n_parents 
#'       The number of parents for each location
#'   * LC
#'       A list defining the precomputed spline structure for each range value
#'       in lc_sequence. Each element of LC is a list with length equal to
#'       length(values). Each sub-element gives the conditional mean of the
#'       spline value as a linear combination "plc" of its parents, and 
#'       the conditional standard deviation as a scalar "elc".
#'   * lc_sequence
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
        correlation_function = function(x1, x2, range) {
            d<- sqrt(sum((x1 - x2)^2))
            if( range == 0 ) return( (d == 0) * 1 )
            (1 + d/range) * exp( -d/range )
        },
        lc_range,
        lc_sequence,
        lc_n = 5,
        graph,
        ...
    ) {
    if( !("matrix" %in% class(x)) ) x<- x |> matrix(ncol = 1)
    n_parents<- min(n_parents, nrow(x))
    
    if( missing(graph) ) {
        graph<- x |> dist() |> distance_matrix_to_dag(k = n_parents)
    }
    if( missing(lc_sequence) & missing(lc_range) ) {
        if( missing(lc_range) ) {
            upper_lc<- find_upper_range(
                x,
                graph,
                correlation_function,
                ...
            )
            lc_range<- c(0, upper_lc)
        }
        lc_sequence<- seq(min(lc_range), max(lc_range), length.out = lc_n)
    }
    lc_sequence<- lc_sequence |> sort()

    spline<- list(
        values = x |> nrow() |> numeric(),
        map = x |> nrow() |> seq(),
        x = x,
        graph = graph,
        n_parents = n_parents,
        LC = list(),
        lc_sequence = 0,
        mixtape = function() {},
        correlation_function = correlation_function
    )
    spline<- spline |> cache_lc(lc_sequence)
    spline
}

mixtape_recorder<- function(nodes) {
    if( requireNamespace("RTMB") ) {
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
#' @param corfun The correlation function
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
find_upper_range<- function(
        x,
        graph,
        corfun,
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
                            S[i, j]<- S[j, i]<- corfun(
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
    target_range<- obj |> 
        lapply(\(x) nlminb(x$par, x$fn)) |>
        lapply(`[[`, "par") |>
        do.call(c, args = _) |>
        quantile(probs = 0.8) |>
        unname()
    return(target_range)
}

#' Precompute spline structure
#' 
#' @param lc_sequence A numeric vector giving the range parameters used to
#'     precompute the spline structure.
#' @param spline A lcspline object.
#' 
#' @return A copy of spline with an updated precomputed spline structure
#' 
#' @export
cache_lc<- function(spline, lc_sequence) {
    spline$lc_sequence<- lc_sequence
    spline$LC<- lc_sequence |>
        lapply(
            compute_single_lc,
            spline$correlation_function,
            spline$graph,
            spline$x
        )
    spline$mixtape<- mixtape_recorder(lc_sequence)
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
#' @param smoothness
#'     A numeric value between 0 and 1 determining the smoothness of the spline.
#' @param height
#'     A positive numeric value proportional to the height of the spline.
#' @param log
#'     If true, returns the log-likelihood.
#' 
#' @return
#'     The log-likelihood of the spline evaluated at x. If x is a simref object,
#'     then additionally the values of x will replaced with simulated values.
#' 
#' @export
dlcspline<- function(
        x, 
        spline, 
        smoothness, 
        height, 
        log = TRUE,
        constrain_smoothness = TRUE
    ) {
    ll<- 0
    if( x |> inherits("simref") ) {
    }
    # Very weak prior to smoothness away from the edges
    ll<- ll + constrain_smoothness * dnorm(qlogis(smoothness), 0, 1.5, TRUE)
    mix<- spline$mixtape(smoothness)
    plc<- spline |>
        _$LC |>
        lapply(`[[`, "plc") |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)

    var<- spline |>
        _$LC |>
        lapply(`[[`, "var") |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)
    if( x |> inherits("simref") ) {
        for( i in spline$graph |> igraph::topo_sort() |> as.numeric() ) {
            p<- spline$graph |> igraph::neighbors(i, "in") |> as.numeric()
            ll<- ll + sum(
                dnorm(
                    x[i],
                    sum(plc[i, p] * x[p]),
                    height * sqrt(var[i]),
                    log = TRUE
                )
            )
        }
    } else {
        mean<- (plc %*% x) |> as.matrix()
        ll<- ll + sum(dnorm(x, mean, height * sqrt(var), log = TRUE))
    }
    if( !log ) ll<- exp(ll)
    return( ll )
}

#' Simulated a lcspline
#' 
#' @param spline
#'     The lcspline to simulate from
#' @param smoothness
#'     A numeric value between 0 and 1 determining the smoothness of the spline.
#' @param height
#'     A positive numeric value proportional to the heigh of the spline.
#' 
#' @return
#'     A copy of the spline with new simulated values.
#' 
#' @export
rlcspline<- function(spline, smoothness, height) {
    mix<- spline$mixtape(smoothness)
    plc<- spline |>
        _$LC |>
        lapply(`[[`, "plc") |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)
    var<- spline |>
        _$LC |>
        lapply(`[[`, "var") |>
        (\(x) Map(`*`, x, as.list(mix)))() |>
        Reduce(`+`, x = _)
    for( i in spline |> _$graph |> igraph::topo_sort() |> as.numeric() ) {
        parents<- spline |> 
            _$graph |> 
            igraph::neighbors(i, "in") |> 
            as.numeric()
        mean<- sum(plc[i, parents] * spline$values[parents])
        spline$values[i]<- rnorm(1, mean, height * sqrt(var[i]))
    }
    return( spline )
}


