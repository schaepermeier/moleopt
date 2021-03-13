#include "mo_descent.h"
#include "utils.h"
#include "vector_utils.h"
#include <cmath>
#include <deque>

/*
 * 
 * Gradient Stuff
 * 
 * 
*/

gradient_fn create_gradient_fn(optim_fn fn,
                               const double_vector& lower,
                               const double_vector& upper,
                               const string& method,
                               double eps_gradient) {
  
  gradient_fn g = [fn, lower, upper, method, eps_gradient](const evaluated_point& point) {
    vector<double_vector> gradients(2);
    int d = point.dec_space.size();
    
    // Initialize gradients
    
    gradients[0] = double_vector(d, 0);
    gradients[1] = double_vector(d, 0);
    
    for (int iter_d = 0; iter_d < d; iter_d++) {
      double_vector lower_d(point.dec_space);
      double_vector upper_d(point.dec_space);
      double_vector fn_lower;
      double_vector fn_upper;
      
      if (method == "twosided") {
        lower_d[iter_d] = lower_d[iter_d] - eps_gradient;
        upper_d[iter_d] = upper_d[iter_d] + eps_gradient;
      } else {
        if (point.dec_space[iter_d] == upper[iter_d]) {
          // point == upper
          upper_d[iter_d] = point.dec_space[iter_d];
          lower_d[iter_d] = point.dec_space[iter_d] - eps_gradient;
        } else if (point.dec_space[iter_d] == lower[iter_d]){
          // point == lower
          upper_d[iter_d] = point.dec_space[iter_d] + eps_gradient;
          lower_d[iter_d] = point.dec_space[iter_d];
        } else {
          // Randomized forward or backward differences
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
  };
  
  return g;
}

double_vector mo_steepest_descent_direction(const vector<double_vector>& gradients) {
  int n_objectives = gradients.size();
  
  if (n_objectives != 2) {
    print("Cannot compute descent direction for n_objectives != 2");
    return double_vector();
  }
  
  // Steepest descent = Highest HV Improvement for given step size
  // Indicated by sqrt(HV)
  
  // double_vector mog = -0.5 * (normalize(gradients[0]) + normalize(gradients[1]))
  //   * sqrt(norm(gradients[0])) * sqrt(norm(gradients[1]));
  
  double_vector norms = {norm(gradients[0]), norm(gradients[1])};
  
  if (norms[0] < 1e-8 ||
      norms[1] < 1e-8 ||
      isnan((sqrt(norms[1]) / sqrt(norms[0]))) ||
      isnan((sqrt(norms[0]) / sqrt(norms[1])))) {
    return 0 * gradients[0];
  }
  
  double_vector mog = -0.5 * (
    (sqrt(norms[1]) / sqrt(norms[0])) * gradients[0] +
    (sqrt(norms[0]) / sqrt(norms[1])) * gradients[1]
  );
  
  return mog;
}

double_vector mo_steepest_descent_direction(const vector<double_vector>& gradients,
                                            const double_vector& ref_point,
                                            const double_vector& current_y) {
  if (dominates(current_y, ref_point)) {
    // Hypervolume gradient direction:
    double_vector imp = ref_point - current_y;

    return -0.5 * (
        sqrt(imp[1]) / sqrt(imp[0]) * gradients[0] +
        sqrt(imp[0]) / sqrt(imp[1]) * gradients[1]
    );
  } else {
    return mo_steepest_descent_direction(gradients);
  }
}

/*
 * 
 * Descent Operators
 * 
 * 
 */

corrector_fn create_two_point_stepsize_descent(const optim_fn& fn,
                                               const gradient_fn& grad_fn,
                                               double eps_initial_step_size,
                                               double eps_descent_direction,
                                               const double_vector& lower,
                                               const double_vector& upper) {

  corrector_fn corr_fn = [fn, grad_fn, eps_initial_step_size, eps_descent_direction, lower, upper]
                         (evaluated_point starting_point, double_vector ref_point, double max_descent) {
    
    // TODO change to parameter
    double max_stepsize = 0.1;
    double scale_factor = 2;
    double armijo_scale = 1e-4;
    
    // Non-monotone search setup
    
    int history_size = 20;
    
    std::deque<double_vector> obj_history;
    
    for (int i = 0; i < history_size; i++) {
      obj_history.push_back(ref_point);
    }
    
    // Setup
    
    evaluated_point current_iterate = starting_point;
    evaluated_point previous_iterate;

    vector<double_vector> gradients = grad_fn(starting_point);
    
    double_vector descent_direction = mo_steepest_descent_direction(gradients);
    descent_direction = project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
    
    double_vector previous_descent_direction = descent_direction;
    
    if (norm(descent_direction) == 0) {
      return starting_point;
    }

    double alpha = eps_initial_step_size / norm(descent_direction);
    
    // First: Some line search to find initial alpha
    
    evaluated_point trial_point = starting_point;
    
    do {
      current_iterate = trial_point;
      
      trial_point.dec_space = ensure_boundary(starting_point.dec_space + alpha * descent_direction,
                                              lower, upper);
      trial_point.obj_space = fn(trial_point.dec_space);
      
      alpha *= scale_factor;
    } while (dominates(trial_point.obj_space, current_iterate.obj_space) &&
             norm(current_iterate.dec_space - starting_point.dec_space) < max_descent &&
             alpha * norm(descent_direction) < max_stepsize);
    
    if (norm(current_iterate.dec_space - starting_point.dec_space) >= max_descent) {
      return current_iterate;
    }
    
    // We increased alpha once too often
    alpha /= scale_factor;
    // Set to starting_point and treat previous step as first iteration
    previous_iterate = starting_point;
    previous_descent_direction = descent_direction;
    
    if (current_iterate.dec_space == starting_point.dec_space) {
      // We did not do any progress, so we can stop here
      return starting_point;
    } else {
      gradients = grad_fn(current_iterate);
      descent_direction = mo_steepest_descent_direction(gradients);
      descent_direction = project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
      int iters = 0;
      
      while (norm(descent_direction) >= eps_descent_direction && alpha > 0 &&
             norm(current_iterate.dec_space - starting_point.dec_space) < max_descent) {
        iters++;
        
        double_vector sk = current_iterate.dec_space - previous_iterate.dec_space;
        double_vector yk = descent_direction - previous_descent_direction;
        
        // Barzilai-Borwein
        double alpha_bb = (dot(sk, sk) / dot(sk, -yk));
        // double alpha_bb = (dot(sk, -yk) / dot(yk, yk));
        double alpha_pos = norm(sk) / norm(yk);

        if (alpha_bb <= 0) {
          alpha = alpha_pos;
        } else {
          alpha = alpha_bb;
        }

        if (alpha == inf || alpha == 0 || isnan(alpha)) {
          break;
        }
        
        alpha = min(alpha, max_stepsize / norm(descent_direction));
        // alpha = max(alpha, eps_initial_step_size / norm(descent_direction));
        
        double_vector expected_improvements = {
          dot(normalize(descent_direction), gradients[0]),
          dot(normalize(descent_direction), gradients[1])
        };
        
        // print(expected_improvements);
        
        // Decreased once too often below
        alpha *= scale_factor;
        
        do {
          alpha /= scale_factor;
          
          trial_point.dec_space = ensure_boundary(current_iterate.dec_space + alpha * descent_direction,
                                        lower, upper);
          trial_point.obj_space = fn(trial_point.dec_space);
        } while (!dominates(trial_point.obj_space - armijo_scale * alpha * expected_improvements + 1e-8, ref_point) &&
                  alpha * norm(descent_direction) >= eps_initial_step_size &&
                  alpha * -expected_improvements[0] >= 1e-8 &&
                  alpha * -expected_improvements[1] >= 1e-8
                   );
        
        // print(alpha);
        
        if (!dominates(trial_point.obj_space - armijo_scale * alpha * expected_improvements + 1e-8, ref_point)) {
          break;
        }
        
        // Update State
        
        previous_iterate = current_iterate;
        previous_descent_direction = descent_direction;
        
        current_iterate = trial_point;
        gradients = grad_fn(current_iterate);
        descent_direction = mo_steepest_descent_direction(gradients);
        descent_direction = project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
        
        // Update ref_point
        
        obj_history.push_back(current_iterate.obj_space);
        obj_history.pop_front();
        
        // ref_point = {0, 0};
        // 
        // for (auto& v : obj_history) {
        //   ref_point = ref_point + v;
        // }
        // 
        // ref_point = ref_point / history_size;
        
        ref_point = {-inf, -inf};

        for (auto& v : obj_history) {
          ref_point = {max(ref_point[0], v[0]), max(ref_point[1], v[1])};
        }

        // ref_point = 0.8 * ref_point + 0.2 * current_iterate.obj_space;

        // print(ref_point);
      }
    
      // print("Iters (MO Descent): " + to_string(iters));
    }
    
    return current_iterate;
  };

  return corr_fn;

}
