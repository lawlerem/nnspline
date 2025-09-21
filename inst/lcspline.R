devtools::load_all()
library(RTMB)
library(sf)
library(Matrix)
library(igraph)
library(tmap)
nx<- 2000
spline_pars<- cbind(
    c(1, 1e-5),
    c(1, 0.1),
    c(1, 0.3),
    c(1, 0.6),
    c(1, 2)
)

geometry<- c(xmin = 0, ymin = 0, xmax = 1, ymax = 1) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_sf(geometry = _) |>
    elhelpers::tesselate(nx) |>
    (\(x) {x$tesselation<- NULL; x})()
spline<- create_nnspline(
        geometry |> sf::st_centroid() |> sf::st_coordinates(),
        graph = geometry |> elhelpers::dagify(adjacency_power = 2),
        covariance_function = function(x1, x2, p) {
            d<- sqrt(sum((x1 - x2)^2))
            variance<- p[[1]]
            range<- p[[2]]
            cov<- variance * (1 + d/range) * exp( -d / range )
            return(cov)
        }
    )

Sigma<- spline_pars |>
    apply(
        2,
        \(x) spline |> update_parameters(x) |> _$covariance_matrix
    )
L<- Sigma |>
    lapply(
        \(S) {
            l<- spline$values |>
                seq_along() |>
                lapply(\(i) c(i, spline$graph |> neighbors(i, "in") |> as.numeric())) |>
                lapply(\(x) S[x, x, drop = FALSE]) |>
                lapply(\(x) {
                    AA<- x[1, 1, drop = FALSE]
                    BA<- x[-1, 1, drop = FALSE]
                    BB<- x[-1, -1, drop = FALSE]
                    AB_BBinv<- solve(BB, BA) |> t()
                    plc<- AB_BBinv |> as.matrix()
                    elc<- (AA - AB_BBinv %*% BA) |> as.numeric() |> sqrt()
                    list(plc = plc, elc = elc)
                })
        }
    )



mix<- c(1, 1, 1, 1, 100) |> (\(x) x / sum(x))()
value<- 0
errors<- rnorm(geometry |> nrow())
for( i in spline$graph |> topo_sort() |> as.numeric() ) {
    parents<- spline$graph |> neighbors(i, "in") |> as.numeric()
    plc<- L |>
        lapply(`[[`, i) |>
        lapply(`[[`, "plc") |>
        Map(`*`, x = _, mix) |>
        Reduce(`+`, x = _)
    elc<- L |>
        lapply(`[[`, i) |>
        lapply(`[[`, "elc") |>
        Map(`*`, x = _, mix) |>
        Reduce(`+`, x = _)
    value[i]<- plc %*% value[parents] + elc * errors[i]
}
geometry$value<- value
tm_shape(geometry) + tm_fill(fill = "value", fill.scale = tm_scale_continuous())