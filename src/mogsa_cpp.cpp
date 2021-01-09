#include "mogsa_cpp.h"
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <cmath>

using namespace std;

double_vector lower;
double_vector upper;

double eps_gradient;
double eps_explore_set;
double eps_initial_step_size;
double max_explore_set;

optim_fn fn;

bool verbose = false;

void log(std::string message) {
  if (verbose) {
    std::cout << message << std::endl;
  }
}

void print(std::string message) {
  std::cout << message << std::endl;
}

void print(double a) {
  std::cout << a << std::endl;
}

void print_vector(double_vector v) {
  for (const auto& el : v) std::cout << el << " ";
  
  std::cout << std::endl;
}

std::default_random_engine generator;
std::normal_distribution<double> distribution;

double random_double() {
  return distribution(generator);
}

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
//     // print_vector(perturbation);
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
    
    double_vector gradient_components = (fn_upper - fn_lower) / length;
    
    gradients[0][iter_d] = gradient_components[0];
    gradients[1][iter_d] = gradient_components[1];
  }
  
  return gradients;
}

double_vector compute_descent_direction(const std::vector<double_vector>& gradients) {
  int n_objectives = gradients.size();
  
  if (n_objectives != 2) {
    print("Cannot compute descent direction for n_objectives != 2");
    return double_vector();
  }
  
  // Local HV Improvement (not nice when objectives very differently scaled):

  // double_vector mog = normalize(gradients[0]) + normalize(gradients[1]);
  // mog = -normalize(mog) * dot(normalize(mog), gradients[0]) * dot(normalize(mog), gradients[1]);
  // 
  // return mog;
  
  // Conventional Definition:
  
  return -0.5 * (normalize(gradients[0]) + normalize(gradients[1]));
}

// std::tuple<evaluated_point, std::string> descend_to_set(evaluated_point current_point) {
// 
//   bool converged = false;
//   int iters = 0;
// 
//   std::string status;
// 
//   double eps_step = eps_initial_step_size;
// 
//   evaluated_point gradient_point = current_point;
//   evaluated_point current_best = current_point;
//   evaluated_point next_point;
// 
//   std::vector<double_vector> current_gradients = compute_gradients(current_point);
//   double_vector descent_direction = compute_descent_direction(current_gradients);
// 
//   while (!converged) {
//     iters++;
// 
//     next_point.dec_space = ensure_boundary(gradient_point.dec_space + eps_step * normalize(descent_direction), lower, upper);
//     next_point.obj_space = fn(next_point.dec_space);
// 
//     if (dominates(next_point.obj_space, current_best.obj_space)) {
//       eps_step *= 2;
// 
//       current_best = next_point;
//     } else {
//       if (gradient_point.dec_space != current_best.dec_space) {
//         gradient_point = current_best;
// 
//         current_gradients = compute_gradients(gradient_point);
//         descent_direction = compute_descent_direction(current_gradients);
//       }
// 
//       if (eps_step <= eps_initial_step_size) {
//         converged = true;
// 
//         // if (dominates(current_best.obj_space, next_point.obj_space)) {
//         if (dominates(current_best.obj_space, next_point.obj_space) ||
//             norm(descent_direction) < 1e-2) {
//           status = "regular";
//         } else {
//           status = "degenerate";
//         }
//       }
// 
//       eps_step /= 2;
//     }
//   }
// 
//   // print(to_string(iters));
// 
//   return {current_best, status};
// }

std::tuple<evaluated_point, std::string> descend_to_set(evaluated_point current_point) {

  int iters = 0;

  std::string status = "regular";

  double eps_step = eps_initial_step_size;

  evaluated_point next_point;

  std::vector<double_vector> next_gradients;
  double_vector next_descent_direction;

  std::vector<double_vector> gradients = compute_gradients(current_point);
  double_vector descent_direction = compute_descent_direction(gradients);

  double smoothing_factor = 0.8;

  while (eps_step >= eps_initial_step_size) {
    iters++;

    next_point.dec_space = ensure_boundary(current_point.dec_space + eps_step * normalize(descent_direction), lower, upper);
    next_point.obj_space = fn(next_point.dec_space);

    next_gradients = compute_gradients(next_point);
    next_descent_direction = compute_descent_direction(next_gradients);

    double factor_eps_step = 1 + dot(next_descent_direction, descent_direction) / (norm(next_descent_direction) * norm(descent_direction));

    if (factor_eps_step < 0.2) {
      factor_eps_step = 0.2;
    }
    
    eps_step *= factor_eps_step;

    for (int i = 0; i < gradients.size(); i++) {
      gradients[i] = smoothing_factor * gradients[i] + (1 - smoothing_factor) * next_gradients[i];
    }

    descent_direction = compute_descent_direction(gradients);
    
    if (!dominates(current_point.obj_space, next_point.obj_space)) {
      current_point.dec_space = next_point.dec_space;
      current_point.obj_space = next_point.obj_space;
    }

  }

  print(to_string(iters));

  return {current_point, status};
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
    
    while (adjusted_explore_set > eps_initial_step_size && dominates(current_point.obj_space, next_point.obj_space)) {
      adjusted_explore_set /= 2;
      next_point.dec_space = current_point.dec_space + adjusted_explore_set * set_direction;
      next_point.dec_space = ensure_boundary(next_point.dec_space, lower, upper);
      next_point.obj_space = fn(next_point.dec_space);
    }
    
    auto [descent_point, descent_status] = descend_to_set(next_point);
    next_point = descent_point;

    if (norm(current_point.dec_space - next_point.dec_space) < 0.1 * eps_explore_set) {
      log("very short step");
    }

    if (next_point.obj_space[objective] < current_point.obj_space[objective] && !strictly_dominates(next_point.obj_space, current_point.obj_space)) {
      // successfully made step in set
      previous_point = current_point;
      current_point = next_point;
      
      trace.push_back(current_point);
    } else {
      if (strictly_dominates(next_point.obj_space, current_point.obj_space) &&
          descent_status != "degenerate") {
        return {next_point, trace};
      } else if (next_point.obj_space[objective] >= current_point.obj_space[objective]) {
        log("SO Optimum reached");
        return {{}, trace};
      } else {
        print("Something else lol");
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

std::tuple<std::vector<std::map<double, evaluated_point>>,
           std::vector<std::tuple<int, int>>> run_mogsa(optim_fn f, std::vector<double_vector> starting_points, double_vector lower_bounds, double_vector upper_bounds,
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size, double maximum_explore_set) {
  eps_gradient = epsilon_gradient;
  eps_explore_set = epsilon_explore_set;
  eps_initial_step_size = epsilon_initial_step_size;
  max_explore_set = maximum_explore_set;
  
  fn = f;
  lower = lower_bounds;
  upper = upper_bounds;
  
  std::vector<std::map<double, evaluated_point>> local_sets;
  std::vector<std::tuple<int, int>> set_transitions;
  
  int starting_points_done = 0;
  
  for (double_vector starting_point : starting_points) {
    starting_points_done++;
    print("Starting point No. " + to_string(starting_points_done));
    
    evaluated_point current_point = {
      starting_point,
      fn(starting_point)
    };

    evaluated_point next_point;
    
    std::vector<std::tuple<evaluated_point, int>> points_to_explore;
    
    auto descent_point = descend_to_set(current_point);
    current_point = descent_point;
    
    points_to_explore.push_back({current_point, -1});

    while(points_to_explore.size() > 0) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      current_point = point_to_explore;
      points_to_explore.pop_back();
      
      // Validate that chosen point does not belong to an already explored set
      
      bool already_explored = false;
      int containing_set = 0;
      
      for (auto& set : local_sets) {
        containing_set++;
        
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

            if (current_point.obj_space[0] >= left_neighbor.obj_space[0] && 
                current_point.obj_space[0] <= right_neighbor.obj_space[0] && 
                current_point.obj_space[1] >= right_neighbor.obj_space[1] && 
                current_point.obj_space[1] <= left_neighbor.obj_space[1]) {
              
              double norm_to_left = norm(current_point.dec_space - left_neighbor.dec_space);
              double norm_to_right = norm(current_point.dec_space - right_neighbor.dec_space);
              double norm_left_right = norm(left_neighbor.dec_space - right_neighbor.dec_space);

              if ((norm_to_left * norm_to_left + norm_to_right * norm_to_right) <= (norm_left_right * norm_left_right + eps_initial_step_size)) {
              // if (norm_to_left < norm_left_right && norm_to_right < norm_left_right) {
                set.insert(std::pair<double, evaluated_point>(current_point.obj_space[0], current_point));
                already_explored = true;
                break;
              }
            }
          }
        }
      }
      
      if (already_explored) {
        set_transitions.push_back({origin_set_id, containing_set});
        print("Skipping: Set already explored");
        continue;
      }
      
      print("Exploring new set");
      print_vector(current_point.dec_space);
      print("Points left: " + to_string(points_to_explore.size()));
      
      std::map<double, evaluated_point> current_set;
      int current_set_id = local_sets.size() + 1;
      
      set_transitions.push_back({origin_set_id, current_set_id});
      
      for (int obj = 0; obj < 2; obj++) {
        const auto [next_point, trace] = explore_efficient_set(current_point, obj);
        
        if (next_point.dec_space.size() != 0) {
          points_to_explore.push_back({next_point, current_set_id});
        }
        
        for (const auto& eval_point : trace) {
          current_set.insert(std::pair<double, evaluated_point>(eval_point.obj_space[0], eval_point));
        }
      }
      
      local_sets.push_back(current_set);
    }
  }
  
  return {local_sets, set_transitions};
}
