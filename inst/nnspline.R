library(nnspline)
library(RTMB)
library(sf)
library(Matrix)
library(igraph)
library(tmap)
nx<- 10000
ny<- 500
w_range<- 0.08
w_sd<- 5
y_sd<- 3


geometry<- c(xmin = 0, ymin = 0, xmax = 1, ymax = 1) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_sf(geometry = _) |>
    elhelpers::tesselate(nx) |>
    (\(x) {x$tesselation<- NULL; x})()
spline<- create_nnspline(
        geometry |> sf::st_centroid() |> sf::st_coordinates(),
        parameters = c(w_sd^2, w_range),
        graph = geometry |> elhelpers::dagify(adjacency_power = 2)
    )
spline<- spline |> rspline()
geometry$value<- spline$values
tm_shape(geometry) + tm_fill(fill = "value", fill.scale = tm_scale_continuous())
y<- geometry |>
    nrow() |>
    sample(ny, replace = TRUE) |>
    (\(idx) {
        data.frame(
            idx = idx,
            y = rnorm(ny, geometry$value[idx], y_sd)
        )
    })()

f<- function(pars) {
    pars |> RTMB::getAll()
    ll<- 0

    spline<- spline |> 
        update_parameters(
            c(log_w_sd, qsmooth) |> exp() |> (\(x) x^c(2, 1))()
        )
    ll<- ll + dspline(w, spline, TRUE)
    ll<- ll + sum(dnorm(y$y, w[y$idx], exp(log_y_sd), TRUE))

    return( -ll )
}
pars<- list(
    qsmooth = -1, 
    log_w_sd = 0,
    w = 0 * spline$values,
    log_y_sd = 0
)
system.time({
    obj<- f |> MakeADFun(pars, random = "w")
    opt<- obj |> with(nlminb(par, fn, gr))
    sdr<- obj |> sdreport(opt$par)
})
