#include <Rcpp.h>
#include "mogsa_cpp.h"
#include "mo_descent.h"
#include "explore_set.h"
#include "utils.h"

using namespace Rcpp;
using namespace std;

optim_fn as_optim_fn(Function fn, long& used_budget, long& max_budget) {
  optim_fn f = [fn, used_budget, max_budget](const double_vector& x) mutable {
    if (max_budget > 0 && used_budget >= max_budget) {
      double_vector retval(2, inf);
      return retval;
    }
    
    NumericVector result = fn(x);
    used_budget++;
    
    // print(used_budget);
    
    return as<double_vector>(result);
  };
  
  return f;
}

corrector_fn as_corrector_fn(Function descent_fn) {
  corrector_fn f = [descent_fn](evaluated_point x, double_vector ref_point, double max_descent) {
    List result = descent_fn(x.dec_space, ref_point);
    
    evaluated_point retval = {
      as<double_vector>(result["dec_space"]),
      as<double_vector>(result["obj_space"])
    };
    
    return retval;
  };
  
  return f;
}

vector<double_vector> rows_to_vectors(NumericMatrix m) {
  vector<double_vector> vectors;
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
  
  vector<double_vector> vectors = rows_to_vectors(m);
  LogicalVector nondominated(n_vectors);
  
  vector<int> idx(n_vectors);
  iota(idx.begin(), idx.end(), 0);
  
  stable_sort(idx.begin(), idx.end(), [&] (int i, int j) {return vectors[i] < vectors [j];});
  
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
List run_mogsa_cpp(
    Function fn,
    NumericMatrix starting_points,
    NumericVector lower_bounds,
    NumericVector upper_bounds,
    double epsilon_gradient,
    double epsilon_explore_set,
    double epsilon_initial_step_size,
    double max_explore_set,
    Nullable<Function> custom_descent_fn = R_NilValue,
    long max_budget = inf,
    std::string logging = "info") {
  
  if (logging == "none") {
    MOLE_LOG_LEVEL = MOLE_LOG_LEVEL_NONE;
  } else if (logging == "debug") {
    MOLE_LOG_LEVEL = MOLE_LOG_LEVEL_DEBUG;
  } else {
    MOLE_LOG_LEVEL = MOLE_LOG_LEVEL_INFO;
  }
  
  /* ========= Setup and run Mogsa ========= */
  
  long used_budget = 0;
  
  optim_fn mo_function = as_optim_fn(fn, used_budget, max_budget);
  double_vector lower = as<double_vector>(lower_bounds);
  double_vector upper = as<double_vector>(upper_bounds);
  
  // Create Gradient of fn
  
  gradient_fn gradient_function = create_gradient_fn(mo_function,
                                           lower,
                                           upper,
                                           "twosided",
                                           epsilon_gradient);
  
  // Create the descent function
  
  corrector_fn descent_function;
  
  if (custom_descent_fn.isNotNull()) {
    Function unpacked_descent_fn(custom_descent_fn);
    descent_function = as_corrector_fn(unpacked_descent_fn);
  } else {
    descent_function = create_two_point_stepsize_descent(mo_function,
                                                         gradient_function,
                                                         1e-8,
                                                         epsilon_initial_step_size,
                                                         0.1,
                                                         2,
                                                         1e-4,
                                                         100,
                                                         1000,
                                                         lower,
                                                         upper);
    
  }
  
  // Create the explore_set function
  
  explore_set_fn explore_set_function = get_explore_set_fn(
    mo_function,
    gradient_function,
    descent_function,
    lower,
    upper,
    epsilon_explore_set,
    max_explore_set,
    45,
    2);
  
  // Run Mogsa
  
  auto [local_sets, set_transitions] = run_mogsa(
            mo_function,
            gradient_function,
            descent_function,
            explore_set_function,
            rows_to_vectors(starting_points),
            lower,
            upper,
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

