#' Create a nnspline object
#'
#' @param x 
#'     A numeric matrix giving the locations (as rows) for which the spline
#'     values are desired
#' @param n_parents 
#'     The number of parents each node should have.
#' @param parameters 
#'     The parameters used in the covariance function.
#' @param noise
#'     A non-negative variance parameter used to model uncorrelated noise.
#' @param covariance_function 
#'     A function used to compute covariances between spline inputs. 
#'     Should have the arguments x1 and x2, vectors giving the spline inputs, 
#'     and p, a vector of parameters.
#' @param graph 
#'     An optional argument. If supplied, it needs to be an igraph object 
#'     representing a directed acyclic graph.
#'
#' @return 
#'     A spline object as a list with the following elements:
#'   * values 
#'       The output values of the spline at x.
#'   * x 
#'       The locations at which the spline is evaluated.
#'   * map 
#'       Mapping indices for the node values. See ?TMB::MakeADFun.
#'   * graph 
#'       An igraph directed acyclic graph represnting the node parent structure.
#'   * n_parents 
#'       The number of parents for each location
#'   * parameters 
#'       The parameters for the covariance function.
#'   * parameter_map 
#'       Mapping indices for the parameters. See ?TMB::MakeADFun.
#'   * noise
#'       A non-negative variance parameter used to model uncorrelated noise.
#'   * covariance_function
#'       The function used to compute covariances between input points.
#'       It should take argument d1 and d2, each of which represents a row
#'       of either spline$x or spline$nodes, and an argument p giving
#'       the parameters of the covariance function.
#'   * covariance_matrix
#'       A partially filled in covariance matrix for the spline nodes. Only the
#'       elements that are needed to either compute the spline likelihood or
#'       update the spline values are filled in, all other entries are set
#'       equal to -Inf.
#'   * pairs
#'       A matrix giving the indices of covariance_matrix that need to be 
#'       filled in.
#' 
#' @export
create_nnspline<- function(
        x,
        n_parents = 4,
        parameters = c(1, 1),
        noise = 0,
        covariance_function = function(x1, x2, p) {
            d<- sqrt(sum((x1 - x2)^2))
            variance<- p[[1]] * p[[2]]^(-3)
            range<- p[[2]]
            cov<- variance * (1 + d / range) * exp( -d / range )
            return(cov)
        },
        graph
    ) {
    if( !("matrix" %in% class(x)) ) x<- x |> matrix(ncol = 1)
    n_parents<- min(n_parents, nrow(x))

    if( missing(graph) ) {
        graph<- x |> 
            dist() |> 
            distance_matrix_to_dag(k = n_parents)
    }
    pairs<- get_pairs(graph)
    covariance_matrix<- matrix(-Inf, nrow = nrow(x), ncol = nrow(x))

    spline<- list(
        values = numeric(nrow(x)),
        x = x,
        map = seq(nrow(x)),
        graph = graph,
        n_parents = n_parents,
        parameters = parameters,
        parameter_map = seq_along(parameters),
        noise = noise,
        covariance_function = covariance_function,
        covariance_matrix = covariance_matrix,
        pairs = pairs
    )
    spline<- spline |> update_spline()
    return(spline)
}

#' Update a spline based on node values and parameters values
#'
#' @param spline 
#'     An nnspline
#' @param ... 
#'     Used to pass arguments onto update_parameters and
#'     update_values.
#'
#' @return 
#'     An updated copy of the spline.
#'
#' @export
update_spline<- function(
        spline,
        parameters,
        values,
        ...
    ) {
    spline<- spline |> 
        update_parameters(...) |>
        update_values(...)

    return(spline)
}

#' @describeIn update_spline 
#'     Update a spline's covariance matrix
#' 
#' @param parameters 
#'     The new parameters for the spline. If missing, will use the pre-existing
#'     parameters stored in the spline.
#' @param noise
#'     The new variance parameter for uncorrelated noise. If missing, will use
#'     the pre-existing noise parameter stored in the spline.
#' 
#' @export
update_parameters<- function(
		spline,
		parameters,
		noise,
		...
	) {
	if( !missing(parameters) ) spline$parameters<- parameters
	if( !missing(noise) ) spline$noise<- noise
	mode<- "numeric"
	if( "advector" %in% class(spline$parameters) && requireNamespace("RTMB") ) {
		mode<- "advector"
		spline$covariance_matrix<- spline$covariance_matrix |> 
            RTMB::AD(force = TRUE)
	}

	covs<- spline$pairs |> 
        apply(
	    	MARGIN = 1,
	    	function(pair) spline$covariance_function(
	    		spline$x[pair[[1]], ],
	    		spline$x[pair[[2]], ],
	    		spline$parameters
	    	),
	    	simplify = FALSE
	    )
	covs<- c |> do.call(covs)
	spline$covariance_matrix[spline$pairs]<- covs
	spline$covariance_matrix[spline$pairs[, c(2, 1)]]<- covs
	diag(spline$covariance_matrix)<- diag(spline$covariance_matrix) + spline$noise

	return(spline)
}

#' @describeIn update_spline 
#'     Update a spline's values
#' 
#' @param values 
#'     The new values for the spline. If missing, will use the pre-existing 
#'     node values stored in the spline.
#' 
#' @export
update_values<- function(
        spline,
        values,
        ...
    ) {
    if( !missing(values) ) spline$values<- values
    mode<- "numeric"
    if( "advector" %in% class(spline$values) && requireNamespace("RTMB") ) {
        mode<- "advector"
        spline$values<- spline$values |> RTMB::AD(force = TRUE)
    }

    return( spline )
}

#' Evaluate the log-likelihood of a spline
#'
#' @param x 
#'     The node values
#' @param spline 
#'     The spline object
#' @param log 
#'     If TRUE, returns the log-likelihood.
#'
#' @return 
#'     The log-likelihood of the spline evaluated at x. If x is a simref
#'     object, then the values of x will be replaced with simulated values.
#'
#' @export
dspline<- function(
        x,
        spline,
        log = TRUE
    ) {
    mode<- "numeric"
    if( "advector" %in% c(class(x), class(spline$covariance_matrix)) && 
            requireNamespace("RTMB") ) mode<- "advector"
    order<- spline$graph |> igraph::topo_sort() |> as.numeric()
    ll<- 0
    for( i in order ) {
        p<- spline$graph |> igraph::neighbors(i, "in") |> as.numeric()
        Sigma<- spline$covariance_matrix[c(i, p), c(i, p), drop = FALSE]
        if( "simref" %in% (x |> class()) ) {
            the_x<- x[c(i, p)]
            dim(the_x)<- c(length(p) + 1, 1)
        } else {
            the_x<- x[c(i, p)] |> cbind()
        }
        cmvn<- conditional_normal(
            the_x,
            (0 * spline$values[c(i, p)]) |> cbind(),
            Sigma
        )
        if( (mode == "advector") | ("simref" %in% (x |> class())) ) {
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
    if( !log ) ll<- ll |> exp()
    return(ll)
}

#' Simulate a spline
#' 
#' @param spline 
#'     The spline to simulate from
#' @param n 
#'     The number of samples
#' 
#' @return 
#'     A list of length n with the simulated splines. If n == 1 then the 
#'     spline will be returned not in a list.
#' 
#' @export
rspline<- function(
        spline,
        n = 1
    ) {
    spline<- spline |> update_parameters()
    values<- matrix(
        0,
        nrow = spline$values |> length(),
        ncol = n
    )
    order<- spline$graph |> igraph::topo_sort() |> as.numeric()
    for( i in order ) {
        p<- spline$graph |> igraph::neighbors(i, "in") |> as.numeric()
        Sigma<- spline$covariance_matrix[c(i, p), c(i, p), drop = FALSE]
        cmvn<- conditional_normal(
            values[c(i, p), , drop = FALSE],
            0 * values[c(i, p), , drop = FALSE],
            Sigma
        )
        values[i, ]<- MASS::mvrnorm(
            n = values |> ncol(),
            mu = 0 * cmvn$mean[, 1],
            Sigma = cmvn$covariance
        )
        values[i, ]<- values[i, ] + cmvn$mean
    }
    sims<- seq(n) |> 
        lapply(
            function(i) spline |> update_values(values[, i])
        )

    if( length(sims) == 1 ) sims<- sims[[1]]
    return(sims)
}
