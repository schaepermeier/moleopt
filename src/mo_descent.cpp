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

gradient_fn create_gradient_fn(const optim_fn& fn,
                               const double_vector& lower,
                               const double_vector& upper,
                               const string& method,
                               double eps_gradient) {
  
  gradient_fn g = [&fn, lower, upper, method, eps_gradient](const evaluated_point& point) mutable {
    vector<double_vector> gradients(2);
    int d = point.dec_space.size();
    
    // Initialize gradients
    
    gradients[0] = double_vector(d, 0);
    gradients[1] = double_vector(d, 0);
    
    double_vector lower_d;
    double_vector upper_d;
    double_vector fn_lower;
    double_vector fn_upper;
    double length;
    double_vector gradient_components;
    
    for (int iter_d = 0; iter_d < d; iter_d++) {
      lower_d = point.dec_space;
      upper_d = point.dec_space;

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
      
      ensure_boundary(lower_d, lower, upper);
      ensure_boundary(upper_d, lower, upper);
      
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
      
      length = norm(lower_d - upper_d);
      
      gradient_components = (fn_upper - fn_lower) / length;
      
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
  
  double_vector norms = {norm(gradients[0]), norm(gradients[1])};
  double factor_grad_1 = sqrt(norms[1] / norms[0]);
  
  if (
      // norms[0] < 1e-6 ||
      // norms[1] < 1e-6 ||
      isnan(factor_grad_1) ||
      isnan(1 / factor_grad_1)) {
    return 0 * gradients[0];
  }
  
  double_vector mog = -0.5 * (
    factor_grad_1 * gradients[0] +
    1 / factor_grad_1 * gradients[1]
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
                                               const double_vector& lower,
                                               const double_vector& upper,
                                               double descent_direction_min,
                                               double descent_step_min,
                                               double descent_step_max,
                                               double descent_scale_factor,
                                               double descent_armijo_factor,
                                               int descent_history_size,
                                               int descent_max_iter) {

  corrector_fn corr_fn = [&fn, &grad_fn, descent_direction_min,
                          descent_step_min, descent_step_max,
                          descent_scale_factor, descent_armijo_factor,
                          descent_history_size, descent_max_iter,
                          lower, upper]
  (const evaluated_point& starting_point, double_vector ref_point, double max_descent) mutable {
                            
    std::deque<double_vector> obj_history;
    obj_history.push_back(ref_point);

    // Setup

    evaluated_point current_iterate = starting_point;
    evaluated_point previous_iterate;

    vector<double_vector> gradients = grad_fn(starting_point);

    double_vector descent_direction = mo_steepest_descent_direction(gradients);
    project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
    double norm_descent_direction = norm(descent_direction);

    double_vector previous_descent_direction = descent_direction;

    if (norm_descent_direction == 0) {
      return starting_point;
    }

    double alpha = descent_step_min / norm_descent_direction;

    // First: Some line search to find initial alpha

    evaluated_point trial_point = starting_point;

    do {
      current_iterate = trial_point;

      trial_point.dec_space = starting_point.dec_space + alpha * descent_direction;
      ensure_boundary(trial_point.dec_space, lower, upper);
      
      trial_point.obj_space = fn(trial_point.dec_space);

      alpha *= descent_scale_factor;
    } while (dominates(trial_point.obj_space, current_iterate.obj_space) &&
             norm(current_iterate.dec_space - starting_point.dec_space) < max_descent &&
             alpha * norm_descent_direction < descent_step_max);

    if (norm(current_iterate.dec_space - starting_point.dec_space) >= max_descent) {
      return current_iterate;
    }

    // We increased alpha once too often
    alpha /= descent_scale_factor;
    // Set to starting_point and treat previous step as first iteration
    previous_iterate = starting_point;
    previous_descent_direction = descent_direction;

    if (current_iterate.dec_space == starting_point.dec_space) {
      // We did not do any progress, so we can stop here
      return starting_point;
    } else {
      gradients = grad_fn(current_iterate);
      descent_direction = mo_steepest_descent_direction(gradients);
      project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
      norm_descent_direction = norm(descent_direction);
      
      int iters = 0;

      while (norm_descent_direction > descent_direction_min && alpha > 0 &&
             norm(current_iterate.dec_space - starting_point.dec_space) < max_descent &&
             iters < descent_max_iter) {
        iters++;

        double_vector sk = current_iterate.dec_space - previous_iterate.dec_space;
        double_vector yk = descent_direction - previous_descent_direction;

        // Barzilai-Borwein
        double alpha_bb = (dot(sk, sk) / dot(sk, -yk));
        double alpha_pos = norm(sk) / norm(yk);

        if (alpha_bb <= 0) {
          alpha = alpha_pos;
        } else {
          alpha = alpha_bb;
        }

        if (alpha == inf || alpha == 0 || isnan(alpha)) {
          break;
        }

        alpha = min(alpha, descent_step_max / norm_descent_direction);
        alpha = max(alpha, descent_step_min / norm_descent_direction);

        double_vector expected_improvements = {
          dot(descent_direction / norm_descent_direction, -gradients[0]),
          dot(descent_direction / norm_descent_direction, -gradients[1])
        };
        
        trial_point.dec_space = current_iterate.dec_space + alpha * descent_direction;
        ensure_boundary(trial_point.dec_space, lower, upper);
        
        trial_point.obj_space = fn(trial_point.dec_space);

        while (!dominates(trial_point.obj_space + descent_armijo_factor * alpha * expected_improvements, ref_point) &&
                alpha > descent_step_min / norm_descent_direction) {
          alpha = max(alpha / descent_scale_factor, descent_step_min / norm_descent_direction);
          
          trial_point.dec_space = current_iterate.dec_space + alpha * descent_direction;
          ensure_boundary(trial_point.dec_space, lower, upper);
          
          trial_point.obj_space = fn(trial_point.dec_space);
        }

        if (!dominates(trial_point.obj_space + descent_armijo_factor * alpha * expected_improvements, ref_point) ||
            (alpha * norm_descent_direction <= descent_step_min && !dominates(trial_point.obj_space, current_iterate.obj_space))) {
          break;
        }

        // print(current_iterate.obj_space - trial_point.obj_space);

        // Update State

        previous_iterate = current_iterate;
        previous_descent_direction = descent_direction;

        current_iterate = trial_point;
        gradients = grad_fn(current_iterate);
        descent_direction = mo_steepest_descent_direction(gradients);
        project_feasible_direction(descent_direction, current_iterate.dec_space, lower, upper);
        norm_descent_direction = norm(descent_direction);

        // Update ref_point

        obj_history.push_back(current_iterate.obj_space);
        if (obj_history.size() > descent_history_size) {
          obj_history.pop_front();
        }

        // ref_point = {0, 0};
        // 
        // for (auto& v : obj_history) {
        //   ref_point = ref_point + v;
        // }
        // 
        // ref_point = ref_point / descent_max_iter;

        ref_point = {-inf, -inf};

        for (auto& v : obj_history) {
          ref_point = {max(ref_point[0], v[0]), max(ref_point[1], v[1])};
        }
        
        // ref_point = 0.8 * ref_point + 0.2 * current_iterate.obj_space;

        // print(ref_point);
      }

      print("Iters (MO Descent): " + to_string(iters));
    }

    return current_iterate;
  };

  return corr_fn;

}
