devtools::load_all()
library(RTMB)
library(sf)
library(Matrix)
library(igraph)
library(tmap)
nx<- 10000
ny<- 2000
lc_sequence<- c(0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
smoothness<- 0.1
height<- 5
mean<- 500
y_sd<- 3

geometry<- c(xmin = 0, ymin = 0, xmax = 1, ymax = 1) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_sf(geometry = _) |>
    elhelpers::tesselate(nx) |>
    (\(x) {x$tesselation<- NULL; x})()
spline<- create_lcspline(
        geometry |> sf::st_centroid() |> sf::st_coordinates(),
        graph = geometry |> elhelpers::dagify(adjacency_power = 2),
        lc_sequence = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
    )
spline<- spline |> rlcspline(smoothness, height)
spline$values |> dlcspline(spline, smoothness, height)
geometry$value<- spline$values
y<- geometry |>
    nrow() |>
    sample(ny, replace = TRUE) |>
    (\(idx) {
        data.frame(
            idx = idx,
            y = rnorm(ny, mean + geometry$value[idx], y_sd)
        )
    })()

f<- function(pars) {
    pars |> RTMB::getAll()
    ll<- 0
    
    ll<- ll + dlcspline(w, spline, plogis(qsmooth), exp(log_height), TRUE)
    ll<- ll + sum(dnorm(y$y, mean + w[y$idx], exp(log_y_sd), TRUE))

    return( -ll )
}
pars<- list(
    qsmooth = 0,
    log_height = 0,
    w = 0 * spline$values,
    mean = 0,
    log_y_sd = 0
)
system.time({
    obj<- f |> MakeADFun(pars, random = "w")
    opt<- obj |> with(nlminb(par, fn, gr))
    sdr<- obj |> sdreport(opt$par)
    sim<- obj$simulate()
})

geometry$estimate<- obj$env$parList(opt$par)$w
geometry$sim<- sim$w
tm_shape(geometry) + 
        tm_fill(
            fill = c("value", "estimate", "sim"), 
            fill.scale = tm_scale_continuous()
        )

