#include <Rcpp.h>
#include "mogsa_cpp.h"

using namespace Rcpp;
using namespace std;

optim_fn as_vector_fn(Function fn) {
  optim_fn f = [fn](double_vector x) {
    NumericVector result = fn(x);
    return as<double_vector>(result);
  };
  
  return f;
}

std::vector<double_vector> rows_to_vectors(NumericMatrix m) {
  std::vector<double_vector> vectors;
  NumericVector row_vector;
  
  for (int i = 0; i < m.nrow(); i++) {
    row_vector = m(i,_);
    
    double_vector v = as<double_vector>(row_vector);
    vectors.push_back(v);
  }
  
  return vectors;
}

// [[Rcpp::export]]
LogicalVector nondominated(NumericMatrix m) {
  int n_vectors = m.nrow();
  
  std::vector<double_vector> vectors = rows_to_vectors(m);
  LogicalVector nondominated(n_vectors);
  
  vector<int> idx(n_vectors);
  iota(idx.begin(), idx.end(), 0);
  
  std::stable_sort(idx.begin(), idx.end(), [&] (int i, int j) {return vectors[i] < vectors [j];});
  
  double best_f2 = vectors[idx[0]][1];
  nondominated[idx[0]] = true;
  
  for (int i = 1; i < n_vectors; i++) {
    if (vectors[idx[i]][1] < best_f2) {
      best_f2 = vectors[idx[i]][1];
      nondominated[idx[i]] = true;
    }
  }
  
  return nondominated;
}

// [[Rcpp::export]]
List run_mogsa_cpp(Function fn, NumericMatrix starting_points, NumericVector lower, NumericVector upper,
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size, double max_explore_set) {
  
  /* ========= Setup and run Mogsa ========= */
  
  auto [local_sets, set_transitions] = run_mogsa(as_vector_fn(fn),
            rows_to_vectors(starting_points),
            as<double_vector>(lower),
            as<double_vector>(upper),
            epsilon_gradient,
            epsilon_explore_set,
            epsilon_initial_step_size,
            max_explore_set);
  
  /* ========= Postprocess the result ========= */
  
  // Convert locally efficient sets to suitable R objects
  
  List sets;
  NumericVector v;
  double d = starting_points.ncol();
  
  for (const auto& set : local_sets) {
    List set_r;
    
    NumericMatrix dec_space(set.size(), d);
    NumericMatrix obj_space(set.size(), (*(set.begin())).second.obj_space.size() );
    int idx = 0;
    
    for (const auto& [key, value] : set) {
      v = wrap(value.dec_space);
      dec_space(idx,_) = v;
      
      v = wrap(value.obj_space);
      obj_space(idx,_) = v;
      
      idx++;
    }
    
    set_r["dec_space"] = dec_space;
    set_r["obj_space"] = obj_space;
    
    sets.insert(sets.size(), set_r);
  }
  
  // Convert set transitions to suitable R objects
  
  NumericMatrix transitions(set_transitions.size(), 2);
  
  int idx = 0;
  
  for (const auto& [from_id, to_id] : set_transitions) {
    transitions(idx, 0) = from_id;
    transitions(idx, 1) = to_id;
    
    idx++;
  }
  
  // Collect return values
  
  List return_values;
  
  return_values["sets"] = sets;
  return_values["transitions"] = transitions;

  return return_values;
}

