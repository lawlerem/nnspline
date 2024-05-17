#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]



// [[Rcpp::export("order_d_matrix")]]
Eigen::VectorXi order_d_matrix(Eigen::MatrixXd &d) {
  Eigen::VectorXi order = Eigen::VectorXi::LinSpaced(
    d.rows(),
    0,
    d.rows() - 1
  );
  int minParent, minChild;
  for(int i = 1; i < d.rows(); i++) {
    // Find least distance
    d.topRightCorner(i, d.cols() - i).minCoeff(&minParent, &minChild);
    minChild += i;
    // Put closest node next
    d.row(i).swap(d.row(minChild));
    d.col(i).swap(d.col(minChild));
    std::iter_swap(order.data() + i, order.data() + minChild);
  }
  return order;
}


struct refSorter {
  refSorter(const Eigen::VectorXd &d) : d_(d) {}
  bool operator () (const int a, const int b) {
    return d_(a) < d_(b);
  }
  const Eigen::VectorXd d_;
};

// [[Rcpp::export("order_adjacency_matrix")]]
Eigen::VectorXi order_adjacency_matrix(Eigen::MatrixXi &m) {
  Eigen::VectorXi order = Eigen::VectorXi::LinSpaced(
    m.rows(),
    0,
    m.rows() - 1
  );
  for(int i = 0; i < m.rows(); i++) {
    Eigen::VectorXi n_parents = m
      .topRightCorner(i, m.cols() - i)
      .colwise()
      .sum();
    int next_vertex;
    n_parents.maxCoeff(&next_vertex);
    next_vertex += i;
    m.row(i).swap(m.row(next_vertex));
    m.col(i).swap(m.col(next_vertex));
    std::iter_swap(order.data() + i, order.data() + next_vertex);
  }
  return order;
}


// [[Rcpp::export("lowest_k")]]
Eigen::VectorXi lowest_k(const Eigen::VectorXd &d, int k) {
  if( k > d.size() ) {
    k = d.size();
  } else {}
  if( k == 0 ) {
    Eigen::VectorXi v(0);
    return v;
  } else {}

  Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(
    d.size(),
    0,
    d.size() - 1
  );
  std::partial_sort(
    ind.data(),
    ind.data() + k,
    ind.data() + ind.size(),
    refSorter(d)
  );
  return ind.segment(0, k);
}


// [[Rcpp::export(".distance_matrix_to_dag")]]
Rcpp::List distance_matrix_to_dag(const Eigen::Map<Eigen::MatrixXd> &d, const int n_neighbours) {
  Eigen::MatrixXd sorted_d = d;
  // Returns order, and sorts ordered_d
  Eigen::VectorXi order = order_d_matrix(sorted_d);

  int dagSize = d.cols();

  Rcpp::List parent_list(dagSize);
  // std::vector<Eigen::VectorXi> edge_list(dagSize);
  parent_list[0] = Eigen::VectorXi(0);
  for(int i = 1; i < parent_list.size(); i++) {
       Eigen::VectorXi p = lowest_k(
        sorted_d.col(i).head(i),
        n_neighbours
      );
      for(int j = 0; j < p.size(); j++) {
        p(j) += 1;
      }
      parent_list[i] = p;
  }

  for(int i = 0; i < order.size(); i++) {
    order(i) += 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("order") = order,
    Rcpp::Named("parent_list") = parent_list
  );
}