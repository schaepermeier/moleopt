#include "mogsa_cpp.h"
#include <iostream>
#include <random>
#include <set>

double_vector lower;
double_vector upper;

double eps_gradient;
double eps_explore_set;
double eps_initial_step_size;

optim_fn fn;

bool verbose = false;

void log(std::string message) {
  if (verbose) {
    std::cout << message << std::endl;
  }
}

std::default_random_engine generator;
std::normal_distribution<double> distribution;

double random_double() {
  return distribution(generator);
}

double inf = std::numeric_limits<double>::infinity();

// std::vector<double_vector> compute_gradients_stochastic(const evaluated_point& point, int n_samples) {
//   std::vector<double_vector> gradients(2);
//   int d = point.dec_space.size();
//   
//   // Initialize gradients
//   
//   // TODO Arbitrary amount of gradients
//   gradients[0] = double_vector(d, 0);
//   gradients[1] = double_vector(d, 0);
//   
//   for (int i = 0; i < n_samples; i++) {
//     // random stochastic perturbation
//     double_vector perturbation = double_vector(d, 0);
//     
//     for (int di = 0; di < d; di++) {
//       perturbation[di] = random_double();
//     }
// 
//     perturbation = normalize(perturbation);
//     
//     // for (const auto& v : perturbation) std::cout << v << " "; std::cout << std::endl;
//     
//     double_vector to_evaluate = ensure_boundary(point.dec_space + eps_gradient * perturbation, lower, upper);
//     double_vector delta_obj = (fn(to_evaluate) - point.obj_space);
//     double_vector delta_dec = (to_evaluate - point.dec_space);
//     
//     gradients[0] = gradients[0] + (delta_obj[0] / eps_gradient * perturbation);
//     gradients[1] = gradients[1] + (delta_obj[1] / eps_gradient * perturbation);
//   }
//   
//   gradients[0] = gradients[0] / n_samples;
//   gradients[1] = gradients[1] / n_samples;
//   
//   return gradients;
// }

std::vector<double_vector> compute_gradients(const evaluated_point& point) {
  std::vector<double_vector> gradients(2);
  int d = point.dec_space.size();
  
  bool twosided = true;
  
  // Initialize gradients

  // TODO Arbitrary amount of gradients
  gradients[0] = double_vector(d, 0);
  gradients[1] = double_vector(d, 0);
  
  for (int iter_d = 0; iter_d < d; iter_d++) {
    double_vector lower_d(point.dec_space);
    double_vector upper_d(point.dec_space);
    double_vector fn_lower;
    double_vector fn_upper;
    
    if (twosided) {
      lower_d[iter_d] = lower_d[iter_d] - eps_gradient;
      upper_d[iter_d] = upper_d[iter_d] + eps_gradient;
    } else {
      if (point.dec_space[iter_d] == upper[iter_d]) {
        // point = upper
        upper_d[iter_d] = point.dec_space[iter_d];
        lower_d[iter_d] = point.dec_space[iter_d] - eps_gradient;
      } else if (point.dec_space[iter_d] == lower[iter_d]){
        // point = lower
        upper_d[iter_d] = point.dec_space[iter_d] + eps_gradient;
        lower_d[iter_d] = point.dec_space[iter_d];
      } else {
        if (random_double() >= 0) {
          upper_d[iter_d] = point.dec_space[iter_d];
          lower_d[iter_d] = point.dec_space[iter_d] - eps_gradient;
        } else {
          upper_d[iter_d] = point.dec_space[iter_d] + eps_gradient;
          lower_d[iter_d] = point.dec_space[iter_d];
        }
      }
    }
    
    lower_d = ensure_boundary(lower_d, lower, upper);
    upper_d = ensure_boundary(upper_d, lower, upper);
    
    if (lower_d == point.dec_space) {
      fn_lower = point.obj_space;
    } else {
      fn_lower = fn(lower_d);
    }
    
    if (upper_d == point.dec_space) {
      fn_upper = point.obj_space;
    } else {
      fn_upper = fn(upper_d);
    }
    
    double length = norm(lower_d - upper_d);
    
    if (length == 0) {
      std::cout << "HELP" << std::endl;
    }
    
    double_vector gradient_components = (fn_upper - fn_lower) / length;
    
    gradients[0][iter_d] = gradient_components[0];
    gradients[1][iter_d] = gradient_components[1];
  }
  
  return gradients;
}

double_vector compute_descent_direction(const std::vector<double_vector>& gradients) {
  // TODO take lower, upper into account!
  int n_objectives = gradients.size();
  
  if (n_objectives != 2) {
    std::cout << "Cannot compute descent direction for n_objectives != 2" << std::endl;
    return double_vector();
  }
  
  return -0.5 * (normalize(gradients[0]) + normalize(gradients[1]));
}

evaluated_point descend_to_set(evaluated_point current_point) {
  
  bool converged = false;
  int iters = 0;
  
  while (!converged) {
    iters++;
    std::vector<double_vector> current_gradients = compute_gradients(current_point);
    double_vector descent_direction = compute_descent_direction(current_gradients);
    
    double eps_step = eps_initial_step_size;
    
    std::pair<double, evaluated_point> a = {0, current_point};
    std::pair<double, evaluated_point> b = {0, current_point};
    std::pair<double, evaluated_point> c = {0, current_point};
    
    evaluated_point new_point;
    
    bool descending_path = true;
    
    while (descending_path) {
      new_point.dec_space = ensure_boundary(current_point.dec_space + eps_step * normalize(descent_direction), lower, upper);
      new_point.obj_space = fn(new_point.dec_space);
      
      a = b;
      b = c;
      c = {eps_step, new_point};
      
      if (dominates(c.second.obj_space, b.second.obj_space)) {
        eps_step *= 2;
      } else {
        descending_path = false;

        if (eps_step <= eps_initial_step_size) {
          converged = true;
        }
      }
    }
    
    if (a.second.dec_space != b.second.dec_space && b.second.dec_space != c.second.dec_space) {
      // Lagrange interpolation (i.e. interpolation w/ quadratic surrogate)
      // Execute until convergence < eps_initial_step_size

      // determine new point for each objective
      // then take the minimum

      std::pair<double, evaluated_point> d;
      double previous_d = 1;

      // while (d.first - previous_d > 1e-6 || d.first - previous_d < -1e-6) {
        previous_d = d.first;
        d.first = inf;

        for (int objective = 0; objective < a.second.obj_space.size(); objective++) {
          double divisor = ((b.first - c.first) * a.second.obj_space[objective] +
                            (c.first - a.first) * b.second.obj_space[objective] +
                            (a.first - b.first) * c.second.obj_space[objective]);
          
          double dividend = ((b.first * b.first - c.first * c.first) * a.second.obj_space[objective] +
                            (c.first * c.first - a.first * a.first) * b.second.obj_space[objective] +
                            (a.first * a.first - b.first * b.first) * c.second.obj_space[objective]);
          
          if (divisor == 0) {
            continue;
          }
          
          double obj_d = 0.5 * dividend / divisor;

          if (obj_d < d.first) {
            d.first = obj_d;
          }
        }

        if (d.first != inf) {
          d.second.dec_space = ensure_boundary(current_point.dec_space + d.first * normalize(descent_direction), lower, upper);
          d.second.obj_space = fn(d.second.dec_space);
          
          // std::cout << "a: " << a.first;
          // std::cout << ", b: " << b.first;
          // std::cout << ", c: " << c.first;
          // std::cout << ", d: " << d.first;
          // std::cout << std::endl;
          
          if (dominates(d.second.obj_space, b.second.obj_space)) {
            b = d;
          }
        }

        // if (a.first < d.first < b.first) {
        //   if (dominates(b.second.obj_space, d.second.obj_space)) {
        //     a = d;
        //   } else {
        //     c = b;
        //     b = d;
        //   }
        // } else {
        //   if (dominates(b.second.obj_space, d.second.obj_space)) {
        //     c = d;
        //   } else {
        //     a = b;
        //     b = d;
        //   }
        // }
      // }

    }
    
    current_point = b.second;
  }
  
  // std::cout << iters << std::endl;

  return current_point;
}

std::tuple<evaluated_point, std::vector<evaluated_point>> explore_efficient_set(evaluated_point current_point, int objective) {
  std::vector<double_vector> current_gradients = compute_gradients(current_point);

  evaluated_point previous_point;
  evaluated_point next_point;
  double_vector set_direction;
  
  std::vector<evaluated_point> trace;
  trace.push_back(current_point);
  
  bool finished = false;
  
  double adjusted_explore_set = eps_explore_set;
  
  while (!finished) {
    if (previous_point.dec_space.size() == 0) {
      // estimate direction using single-objective gradient
      // only make a very small step in this case
      set_direction = 0.1 * -normalize(current_gradients[objective]);
    } else {
      set_direction = normalize(current_point.dec_space - previous_point.dec_space);
    }
    
    adjusted_explore_set *= 2;
    if (adjusted_explore_set > eps_explore_set) {
      adjusted_explore_set = eps_explore_set;
    }
    
    next_point.dec_space = current_point.dec_space + adjusted_explore_set * set_direction;
    next_point.dec_space = ensure_boundary(next_point.dec_space, lower, upper);
    next_point.obj_space = fn(next_point.dec_space);
    
    while (adjusted_explore_set > eps_gradient && dominates(current_point.obj_space, next_point.obj_space)) {
      adjusted_explore_set /= 2;
      next_point.dec_space = current_point.dec_space + adjusted_explore_set * set_direction;
      next_point.dec_space = ensure_boundary(next_point.dec_space, lower, upper);
      next_point.obj_space = fn(next_point.dec_space);
    }
    
    next_point = descend_to_set(next_point);
    
    if (norm(current_point.dec_space - next_point.dec_space) < 0.1 * eps_explore_set) {
      log("very short step");
    }

    if (next_point.obj_space[objective] < current_point.obj_space[objective] && !strictly_dominates(next_point.obj_space, current_point.obj_space)) {
      // successfully made step in set
      previous_point = current_point;
      current_point = next_point;
      
      trace.push_back(current_point);
    } else {
      if (strictly_dominates(next_point.obj_space, current_point.obj_space)) {
        return {next_point, trace};
      } else if (next_point.obj_space[objective] >= current_point.obj_space[objective]) {
        log("SO Optimum reached");
        return {{}, trace};
      } else {
        return {{}, trace};
      }
      // either hit single-objective optimum, gone over ridge or some error
      finished = true;
    }
  }
  
  // empty vector if finished
  // current point if gone over ridge
  return {{}, trace};
}

std::vector<std::map<double, evaluated_point>> run_mogsa(optim_fn f, std::vector<double_vector> starting_points, double_vector lower_bounds, double_vector upper_bounds,
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size) {
  eps_gradient = epsilon_gradient;
  eps_explore_set = epsilon_explore_set;
  eps_initial_step_size = epsilon_initial_step_size;
  fn = f;
  lower = lower_bounds;
  upper = upper_bounds;
  
  std::vector<std::map<double, evaluated_point>> local_sets;
  
  for (double_vector starting_point : starting_points) {
    evaluated_point current_point = {
      starting_point,
      fn(starting_point)
    };
    
    // std::vector<double_vector> current_gradients_fd = compute_gradients(current_point);
    // for (const auto& v : current_gradients_fd[0]) std::cout << v << " "; std::cout << std::endl;
    // 
    // std::vector<double_vector> current_gradients_stoch = compute_gradients_stochastic(current_point, 20);
    // for (const auto& v : current_gradients_stoch[0]) std::cout << v << " "; std::cout << std::endl;
    
    evaluated_point next_point;
    
    std::vector<evaluated_point> points_to_explore;
    
    current_point = descend_to_set(current_point);
    
    points_to_explore.push_back(current_point);
    
    while(points_to_explore.size() > 0) {
      current_point = points_to_explore.back();
      points_to_explore.pop_back();
      
      // Validate that chosen point does not belong to an already explored set
      
      bool already_explored = false;
      
      for (auto& set : local_sets) {
        double lower_f1 = (*(set.begin())).second.obj_space[0];
        double upper_f2 = (*(set.begin())).second.obj_space[1];
        double upper_f1 = (*(set.rbegin())).second.obj_space[0];
        double lower_f2 = (*(set.rbegin())).second.obj_space[1];
        
        if (current_point.obj_space[0] >= lower_f1 && 
            current_point.obj_space[0] <= upper_f1 && 
            current_point.obj_space[1] >= lower_f2 && 
            current_point.obj_space[1] <= upper_f2) {
          
          std::map<double, evaluated_point>::iterator it_lower = set.lower_bound(current_point.obj_space[0]);
  
          if ((it_lower != set.end()) && (it_lower != set.begin())) {
            auto& right_neighbor = (*it_lower).second;
            auto& left_neighbor = (*(--it_lower)).second;
            
            // std::cout << current_point.obj_space[0] << " " << current_point.obj_space[1] << std::endl;
            // std::cout << left_neighbor.obj_space[0] << " " << left_neighbor.obj_space[1] << std::endl;
            // std::cout << right_neighbor.obj_space[0] << " " << right_neighbor.obj_space[1] << std::endl;
  
            if (current_point.obj_space[0] >= left_neighbor.obj_space[0] && 
                current_point.obj_space[0] <= right_neighbor.obj_space[0] && 
                current_point.obj_space[1] >= right_neighbor.obj_space[1] && 
                current_point.obj_space[1] <= left_neighbor.obj_space[1]) {
              
              double norm_to_left = norm(current_point.dec_space - left_neighbor.dec_space);
              double norm_to_right = norm(current_point.dec_space - right_neighbor.dec_space);
              double norm_left_right = norm(left_neighbor.dec_space - right_neighbor.dec_space);
              
              // std::cout << norm_to_left << std::endl;
              // std::cout << norm_to_right << std::endl;
              // std::cout << norm_left_right << std::endl;
              
              // if ((norm_to_left * norm_to_left + norm_to_right * norm_to_right) < (norm_left_right * norm_left_right)) {
              if (norm_to_left < norm_left_right && norm_to_right < norm_left_right) {
                set.insert(std::pair<double, evaluated_point>(current_point.obj_space[0], current_point));
                already_explored = true;
                break;
              }
            }
          }
        }
      }
      
      if (already_explored) {
        std::cout << "Skipping point, set already explored" << std::endl;
        continue;
      }
      
      std::cout << "Exploring new set" << std::endl;
      for (const auto& v : current_point.dec_space) std::cout << v << " "; std::cout << std::endl;
      std::cout << "Points left: " << points_to_explore.size() << std::endl;
      
      std::map<double, evaluated_point> current_set;
      
      for (int obj = 0; obj < 2; obj++) {
        const auto [next_point, trace] = explore_efficient_set(current_point, obj);
        
        if (next_point.dec_space.size() != 0) {
          points_to_explore.push_back(next_point);
        }
        
        for (const auto& eval_point : trace) {
          current_set.insert(std::pair<double, evaluated_point>(eval_point.obj_space[0], eval_point));
        }
      }
      
      local_sets.push_back(current_set);
    }
  }
  
  return local_sets;
}
