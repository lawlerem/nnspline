#' Create a nnspline object
#'
#' @param x 
#'     A numeric matrix giving the locations (as rows) for which the spline
#'     values are desired
#' @param nodes 
#'     A numeric matrix giving the spline knot locations as rows.
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
#' @param node_graph 
#'     An optional argument. If supplied, it needs to be an igraph object 
#'     representing a directed acyclic graph.
#'
#' @return 
#'     A multispline object as a list with the following elements:
#'   * values 
#'       The output values of the spline at x.
#'   * x 
#'       The locations at which the spline is evaluated.
#'   * x_parents 
#'       The nodes which are parents of each x. If a parent is a negative value
#'       then the node value of that parent is taken to be the corresponding
#'       x value.
#'   * node_values 
#'       The output values of the spline at the nodes.
#'   * node_values_map 
#'       Mapping indices for the node values. See ?TMB::MakeADFun.
#'   * nodes 
#'       The locations of the spline nodes.
#'   * node_graph 
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
#'   * x_covariance 
#'       A list giving the covariance between x its parents. The first element
#'       of each is the covariance of x with itself.
#'   * node_covariance 
#'       A partially filled in covariance matrix for the spline nodes. Only the
#'       elements that are needed to either compute the spline likelihood or
#'       update the spline values are filled in, all other entries are set
#'       equal to -Inf.
#'   * node_pairs
#'       A matrix giving the indices of node_covariance that need to be 
#'       filled in.
#' 
#' @export
create_nnspline<- function(
		x,
		nodes = x,
		n_parents = 4,
		parameters = c(1, 1),
		noise = 0,
		covariance_function = function(x1, x2, p) {
			d<- sqrt(sum((x1 - x2)^2))
			variance<- p[[1]]
			range<- p[[2]]
			cov<- variance * (1 + d / range) * exp( -d / range )
			return(cov)
		},
		node_graph
	) {
	if( !("matrix" %in% class(x)) ) x<- x |> matrix(ncol = 1)
	if( !("matrix" %in% class(nodes)) ) nodes<- nodes |> matrix(ncol = 1)
	n_parents<- min(n_parents, nrow(nodes))

	t_nodes<- nodes |> t()
	x_parents<- x |> 
        apply(
	    	MARGIN = 1,
	    	function(row) {
	    		d<- t_nodes - row
	    		d<- d^2 |> colSums() |> sqrt()
	    		d_sort<- d |> sort(partial = n_parents)
	    		p<- (d <= d_sort[n_parents]) |> which()
	    		# The negative sign means that x is EQUAL to node[p,]
	    		if( (d[p] <= .Machine$double.eps) |> any() ) {
                    p<- -p[d[p] |> which.min()]
                }
	    		return(p)
	    	},
	    	simplify = FALSE
	    )
	x_covariance<- x_parents |>
        lapply(
	    	function(p) {
	    		if( (length(p) == 1) && (p < 0) ) return(NA)
	    		return(numeric(length(p)))
	    	}
	    )

	if( missing(node_graph) ) {
		node_graph<- nodes |> 
            dist() |> 
            distance_matrix_to_dag(k = n_parents)
	}
	node_pairs<- get_pairs(node_graph, x_parents)
	node_covariance<- matrix(-Inf, nrow = nrow(nodes), ncol = nrow(nodes))

	spline<- list(
		values = numeric(nrow(x)),
		x = x,
		x_parents = x_parents,
		node_values = numeric(nrow(nodes)),
		node_values_map = seq(nrow(nodes)),
		nodes = nodes,
		node_graph = node_graph,
		n_parents = n_parents,
		parameters = parameters,
		parameter_map = seq_along(parameters),
		noise = noise,
		covariance_function = covariance_function,
		x_covariance = x_covariance,
		node_covariance = node_covariance,
		node_pairs = node_pairs
	)
	spline<- spline |> update_spline()
	return(spline)
}

#' Update a spline based on node values and parameters values
#'
#' @param spline 
#'     An nnspline
#' @param ... 
#'     Used to pass arguments onto update_spline_covariance and
#'     update_spline_values.
#'
#' @return 
#'     An updated copy of the spline.
#'
#' @export
update_spline<- function(
		spline,
		...
	) {
	spline<- spline |> 
        update_spline_covariance(...) |>
        update_spline_values(...)

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
#' @param only_node_covariance 
#'     Logical. If TRUE, the cross-covariance matrix between x and nodes will 
#'     not be updated.
#' 
#' @export
update_spline_covariance<- function(
		spline,
		parameters,
		noise,
		only_node_covariance = FALSE,
		...
	) {
	if( !missing(parameters) ) spline$parameters<- parameters
	if( !missing(noise) ) spline$noise<- noise
	mode<- "numeric"
	if( "advector" %in% class(spline$parameters) && requireNamespace("RTMB") ) {
		mode<- "advector"
		spline$x_covariance<- spline$x_covariance |> 
            lapply(
		    	RTMB::AD,
		    	force = TRUE
		    )
		spline$node_covariance<- spline$node_covariance |> 
            RTMB::AD(force = TRUE)
	}

	covs<- spline$node_pairs |> 
        apply(
	    	MARGIN = 1,
	    	function(pair) spline$covariance_function(
	    		spline$nodes[pair[[1]], ],
	    		spline$nodes[pair[[2]], ],
	    		spline$parameters
	    	),
	    	simplify = FALSE
	    )
	covs<- c |> do.call(covs)
	spline$node_covariance[spline$node_pairs]<- covs
	spline$node_covariance[spline$node_pairs[, c(2, 1)]]<- covs
	diag(spline$node_covariance)<- diag(spline$node_covariance) + spline$noise
	if( only_node_covariance ) return(spline)


	spline$x_covariance<- seq(length(spline$x_parents)) |>
        lapply(
	    	function(xi) {
	    		p<- spline$x_parents[[xi]]
	    		if( (length(p) == 1) && (p < 0) ) return(NA)
	    		self_cov<- spline$covariance_function(
	    			spline$x[xi, ],
	    			spline$x[xi, ],
	    			spline$parameters
	    		) + spline$noise
	    		parent_cov<- p |> 
                    sapply(
	    		    	function(nodei) spline$covariance_function(
	    		    		spline$x[xi, ],
	    		    		spline$nodes[nodei, ],
	    		    		spline$parameters
	    		    	),
	    		    	simplify = FALSE
	    		    )
	    		parent_cov<- c |> do.call(parent_cov)
	    		return(c(self_cov, parent_cov))
	    	}
	    )

	return(spline)
}

#' @describeIn update_spline 
#'     Update a spline's values
#' 
#' @param node_values 
#'     The new new values for the spline. If missing, will use the pre-existing 
#'     node values stored in the spline.
#' 
#' @export
update_spline_values<- function(
		spline,
		node_values,
		...
	) {
	if( !missing(node_values) ) spline$node_values<- node_values
	mode<- "numeric"
	if( "advector" %in% class(spline$node_values) && requireNamespace("RTMB") ) {
		mode<- "advector"
		spline$values<- spline$values |> RTMB::AD(force = TRUE)
	}

	ans<- seq(nrow(spline$x)) |> 
        sapply(
	    	function(xi) {
	    		p<- spline$x_parents[[xi]]
	    		if( (length(p) == 1) && (p < 0) ) {
                    return( spline$node_values[abs(p)] )
                }

	    		Sigma<- spline$node_covariance[
	    			c(1, p),
	    			c(1, p)
	    		]
	    		Sigma[1, ]<- spline$x_covariance[[xi]]
	    		Sigma[, 1]<- spline$x_covariance[[xi]]
	    		cmvn<- conditional_normal(
	    			x = c(spline$values[xi], spline$node_values[p]) |> cbind(),
	    			joint_mean = numeric(length(p) + 1) |> cbind(),
	    			joint_covariance = Sigma,
	    			predicted = 1,
	    			predictors = seq_along(p) + 1
	    		)
	    		return( cmvn$mean )
	    	},
	    	simplify = FALSE
	    )
	spline$values<- c |> do.call(ans)
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
	if( "advector" %in% c(class(x), class(spline$node_covariance)) && 
            requireNamespace("RTMB") ) mode<- "advector"
	node_order<- spline$node_graph |> igraph::topo_sort() |> as.numeric()
	ll<- 0
	for( i in node_order ) {
		p<- spline$node_graph |> igraph::neighbors(i, "in") |> as.numeric()
		Sigma<- spline$node_covariance[c(i, p), c(i, p), drop = FALSE]
		if( "simref" %in% (x |> class()) ) {
			the_x<- x[c(i, p)]
			dim(the_x)<- c(length(p) + 1, 1)
		} else {
			the_x<- x[c(i, p)] |> cbind()
		}
		cmvn<- conditional_normal(
			the_x,
			(0 * spline$node_values[c(i, p)]) |> cbind(),
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
	spline<- spline |> update_spline_covariance()
	node_values<- matrix(
		0,
		nrow = spline$node_values |> length(),
		ncol = n
	)
	node_order<- spline$node_graph |> igraph::topo_sort() |> as.numeric()
	for( i in node_order ) {
		p<- spline$node_graph |> igraph::neighbors(i, "in") |> as.numeric()
		Sigma<- spline$node_covariance[c(i, p), c(i, p), drop = FALSE]
		cmvn<- conditional_normal(
			node_values[c(i, p), , drop = FALSE],
			0 * node_values[c(i, p), , drop = FALSE],
			Sigma
		)
		node_values[i, ]<- MASS::mvrnorm(
			n = node_values |> ncol(),
			mu = 0 * cmvn$mean[, 1],
			Sigma = cmvn$covariance
		)
		node_values[i, ]<- node_values[i, ] + cmvn$mean
	}
	sims<- seq(n) |> 
        lapply(
	    	function(i) spline |> update_spline_values(node_values[, i])
	    )

	if( length(sims) == 1 ) sims<- sims[[1]]
	return(sims)
}


#' Evaluate a nnspline
#' 
#' @param x 
#'     The input values. If x is a vector it is treated as a column vector. The 
#'     dimensions of x should be equal to c(DIM, ncol(spline$x)). If 
#'     ncol(spline$x) == 1, then the dimensions can drop the last dimension and 
#'     be just DIM.
#' @param spline 
#'     The spline
#' @param index 
#'     Logical. If TRUE, returns the index instead of the value.
#' 
#' @return 
#'     The value(s) of the spline at x, or its indices if index is TRUE. 
#'     The dimension of the return is equal to DIM.
#'
#' @export
nns<- function(
		x,
		spline,
		index = FALSE
	) {
	# If x is a vector, make it a matrix with 1 column
	if( x |> dim() |> is.null() ) x<- x |> array(dim = c(x |> length(), 1))
	if( (x |> dim() |> tail(1)) != (spline$x |> ncol()) ) {
		if( (spline$x |> ncol()) == 1 ) {
			x<- x |> array(dim = c(x |> dim(), 1))
		} else {
			stop("The last dimension of x must be equal to the number of columns of spline$x.")
		}
	}
	spline_x<- spline$x |> t()
	idx<- x |>
        apply(
	    	MARGIN = seq(length(dim(x)) - 1),
	    	function(row) {
	    			if( row |> is.na() |> any() ) return( NA )
	    			d<- spline_x - row 
	    			d<- d^2 |> colSums() |> sqrt()
	    			idx<- d |> which.min()
	    			return( idx )
	    	}
	    )
	
	ans<- idx |> c()
	if( !index ) {
		if( "advector" %in% (spline$parameters |> class()) ) {
			ans<- idx |> RTMB::AD()
			ans[idx |> is.na()]<- NA |> RTMB::AD()
			ans[!(idx |> is.na())]<- spline$values[idx[!(idx |> is.na())]]
		} else {
			ans<- spline$values[idx]
		}
	}
	ans<- ans |> array(dim = x |> dim() |> head(-1))
	return( ans )
}

