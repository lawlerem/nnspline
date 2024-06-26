# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

order_d_matrix <- function(d) {
    .Call(`_nnspline_order_d_matrix`, d)
}

order_adjacency_matrix <- function(m) {
    .Call(`_nnspline_order_adjacency_matrix`, m)
}

lowest_k <- function(d, k) {
    .Call(`_nnspline_lowest_k`, d, k)
}

.distance_matrix_to_dag <- function(d, n_neighbours) {
    .Call(`_nnspline_distance_matrix_to_dag`, d, n_neighbours)
}

