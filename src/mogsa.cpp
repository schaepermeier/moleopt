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
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size) {
  
  auto local_sets = run_mogsa(as_vector_fn(fn),
            rows_to_vectors(starting_points),
            as<double_vector>(lower),
            as<double_vector>(upper),
            epsilon_gradient,
            epsilon_explore_set,
            epsilon_initial_step_size);
  
  List return_values;
  List sets;

  double d = starting_points.ncol();
  
  NumericVector v;
  
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
  
  return_values["sets"] = sets;

  return return_values;
}

// double vector_length(NumericVector vector) {
//   return(sqrt(sum(vector * vector)));
// }
// 
// bool dominates(NumericVector a, NumericVector b) {
//   return(is_true(all(a <= b)) && is_true(any(a < b)));
// }
// 
// NumericVector ensure_boundary(NumericVector vector, NumericVector lower, NumericVector upper) {
//   return pmin(pmax(vector, lower), upper);
// }
// 
// NumericVector normalize(NumericVector vector) {
//   double length = vector_length(vector);
//   
//   if (length == 0) {
//     return vector;
//   } else {
//     return vector / length;
//   }
// }
// 
// // [[Rcpp::export]]
// std::vector<NumericVector> compute_gradients(Function fn, NumericVector point, NumericVector lower, NumericVector upper, double epsilon_gradient) {
//   std::vector<NumericVector> gradients(2);
//   int d = point.size();
//   
//   // Initialize gradients
//   NumericVector zero(d, 0.0);
//   
//   // TODO Arbitrary amount of gradients
//   gradients[0] = NumericVector::import(zero.begin(), zero.end());
//   gradients[1] = NumericVector::import(zero.begin(), zero.end());
//   
//   for (int iter_d = 0; iter_d < d; iter_d++) {
//     NumericVector lower_d(point.begin(), point.end());
//     NumericVector upper_d(point.begin(), point.end());
//     
//     // TODO Ensure in-bounds
//     lower_d(iter_d) = lower_d(iter_d) - epsilon_gradient;
//     upper_d(iter_d) = upper_d(iter_d) + epsilon_gradient;
//     
//     lower_d = ensure_boundary(lower_d, lower, upper);
//     upper_d = ensure_boundary(upper_d, lower, upper);
//     
//     double length = vector_length(lower_d - upper_d);
//     
//     NumericVector fn_lower = as<NumericVector>(fn(lower_d));
//     NumericVector fn_upper = as<NumericVector>(fn(upper_d));
//     
//     NumericVector gradient_components = (fn_upper - fn_lower) / length;
//     
//     gradients[0](iter_d) = gradient_components[0];
//     gradients[1](iter_d) = gradient_components[1];
//   }
//   
//   return gradients;
// }
// 
// // [[Rcpp::export]]
// NumericVector compute_descent_direction(std::vector<NumericVector> gradients) {
//   // TODO take lower, upper into account!
//   int n_objectives = gradients.size();
//   
//   if (n_objectives != 2) {
//     std::cout << "Cannot compute descent direction for n_objectives != 2";
//     return NumericVector::create();
//   }
//   
//   return -0.5 * (normalize(gradients[0]) + normalize(gradients[1]));
// }
// 
// // [[Rcpp::export]]
// NumericVector descend_to_set(Function fn, NumericVector current_point, NumericVector lower, NumericVector upper,
//                              double epsilon_gradient, double epsilon_initial_step_size) {
//   double epsilon_step = epsilon_initial_step_size;
//   
//   std::vector<NumericVector> current_gradients = compute_gradients(fn, current_point, lower, upper, epsilon_gradient);
//   NumericVector current_fn = fn(current_point);
//   NumericVector descent_direction = compute_descent_direction(current_gradients);
//   
//   int successes = 0;
//   
//   while (epsilon_step >= epsilon_initial_step_size) {
//     NumericVector next_point = current_point + epsilon_step * normalize(descent_direction);
//     next_point = ensure_boundary(next_point, lower, upper);
//     NumericVector next_fn = fn(next_point);
//     
//     if (dominates(next_fn, current_fn)) {
//       current_point = next_point;
//       current_fn = next_fn;
//       epsilon_step *= 2;
//       
//       successes++;
//     } else {
//       if (successes == 0) {
//         epsilon_step /= 2;
//       } else {
//         current_gradients = compute_gradients(fn, current_point, lower, upper, epsilon_gradient);
//         descent_direction = compute_descent_direction(current_gradients);
//         epsilon_step = epsilon_initial_step_size;
//         
//         successes = 0;
//       }
//     }
//   }
//   
//   return current_point;
// }
// 
// NumericVector explore_efficient_set(Function fn, NumericVector current_point, int objective, NumericVector lower, NumericVector upper,
//                                     double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size) {
//   std::vector<NumericVector> current_gradients = compute_gradients(fn, current_point, lower, upper, epsilon_gradient);
//   NumericVector current_fn = fn(current_point);
//   
//   NumericVector previous_point;
//   NumericVector next_point;
//   NumericVector next_fn;
//   NumericVector set_direction;
//   bool finished = false;
//   
//   while (!finished) {
//     if (previous_point.size() == 0) {
//       // estimate direction using single-objective gradient
//       set_direction = -normalize(current_gradients[objective]);
//     } else {
//       set_direction = normalize(current_point - previous_point);
//     }
//     
//     next_point = current_point + epsilon_explore_set * set_direction;
//     next_point = ensure_boundary(next_point, lower, upper);
//     next_point = descend_to_set(fn, next_point, lower, upper, epsilon_gradient, epsilon_initial_step_size);
//     
//     next_fn = fn(next_point);
//     
//     if (next_fn[objective] < current_fn[objective] && !dominates(next_fn, current_fn)) {
//       // successfully made step in set
//       previous_point.assign(current_point.begin(), current_point.end());
//       current_point.assign(next_point.begin(), next_point.end());
//       current_fn.assign(next_fn.begin(), next_fn.end());
//     } else {
//       if (dominates(next_fn, current_fn)) {
//         return next_point;
//       } else {
//         return NumericVector::create();
//       }
//       // either hit single-objective optimum gone over ridge or some error
//       finished = true;
//     }
//   }
//   
//   // empty vector if finished
//   // current point if gone over ridge
//   return current_point;
// }
// 
// // [[Rcpp::export]]
// List run_mogsa(Function fn, NumericVector starting_point, NumericVector lower, NumericVector upper,
//                double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size) {
//   // results["test"] = 0.0;
//   List results;
// 
//   NumericVector current_point(starting_point);
//   NumericVector next_point;
//   
//   std::vector<NumericVector> points_to_explore;
// 
//   current_point = descend_to_set(fn, current_point, lower, upper, epsilon_gradient, epsilon_initial_step_size);
//   
//   points_to_explore.push_back(current_point);
//   
//   while(points_to_explore.size() > 0) {
//     current_point = points_to_explore.back();
//     points_to_explore.pop_back();
//     
//     std::cout << "Exploring new set" << std::endl;
//     print(current_point);
//     std::cout << "Points left: " << points_to_explore.size() << std::endl;
//     
//     for (int obj = 0; obj < 2; obj++) {
//       next_point = explore_efficient_set(fn, current_point, obj, lower, upper, epsilon_gradient, epsilon_explore_set, epsilon_initial_step_size);
//       
//       if (next_point.size() != 0) {
//         points_to_explore.push_back(next_point);
//       }
//     }
//   }
//   
//   return results;
// }
