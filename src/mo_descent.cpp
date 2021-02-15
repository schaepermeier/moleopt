#include "mo_descent.h"
#include "utils.h"
#include "vector_utils.h"
#include <cmath>

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
  
  double_vector mog = -0.5 * (normalize(gradients[0]) + normalize(gradients[1]))
    * sqrt(norm(gradients[0])) * sqrt(norm(gradients[1]));
  
  return mog;
}

double_vector mo_steepest_descent_direction(const vector<double_vector>& gradients,
                                            const double_vector& ref_point,
                                            const double_vector& current_y) {
  double current_improvement = sqrt(max(0.0, ref_point[0] - current_y[0])) *
                               sqrt(max(0.0, ref_point[1] - current_y[1]));
  
  double_vector mog;
  
  if (current_improvement > 0) {
    // Hypervolume gradient direction:
    
    mog = (ref_point[1] - current_y[1]) * gradients[0] +
          (ref_point[0] - current_y[0]) * gradients[1];
    
    // Normalize and get descent direction
    mog = -mog / (2 * current_improvement);
  } else {
    mog = mo_steepest_descent_direction(gradients);
  }
  
  return mog;
}

/*
 * 
 * Hypervolume Stuff
 * 
 * 
 */

optim_fn create_improvement_function(const optim_fn& fn, double_vector ref_point) {
  if (ref_point.size() != 2) {
    throw "Improvement only supported for dim == 2";
  }
  
  optim_fn hv_fn = [fn, ref_point](double_vector dec_space) {
    double_vector val = {
      compute_improvement(ref_point, fn(dec_space))
    };
    
    return val;
  };
  
  return hv_fn;
}

/*
 * 
 * Descent Operators
 * 
 * 
 */

corrector_fn create_armijo_descent_corrector(const optim_fn& fn,
                                             const gradient_fn& grad_fn,
                                             double eps_initial_step_size,
                                             const double_vector& lower,
                                             const double_vector& upper) {
  
  corrector_fn corr_fn = [fn, grad_fn, eps_initial_step_size, lower, upper](evaluated_point current_point, double_vector ref_point) {
    evaluated_point trial_point;
    
    double_vector descent_direction = {1}; // some default value for the norm to be non-zero
    vector<double_vector> gradients;
    
    double alpha = 2 * eps_initial_step_size;
    double factor = 2;
    int iters = 0;
    
    while (alpha > eps_initial_step_size) {
      iters++;
      
      alpha = eps_initial_step_size;
      
      double current_improvement = compute_improvement(current_point.obj_space, ref_point);
      
      double next_point_improvement = current_improvement;
      evaluated_point next_point = current_point;
      
      gradients = grad_fn(current_point);
      
      // HV gradient direction
      descent_direction = mo_steepest_descent_direction(gradients, ref_point, current_point.obj_space);
      
      bool improving = true;
      
      while (improving) {
        trial_point.dec_space = ensure_boundary(current_point.dec_space + alpha * normalize(descent_direction), lower, upper);
        trial_point.obj_space = fn(trial_point.dec_space);
        
        double trial_point_improvement = compute_improvement(trial_point.obj_space, ref_point);
        
        if (trial_point_improvement > next_point_improvement) {
          next_point = trial_point;
          next_point_improvement = trial_point_improvement;
          
          alpha *= factor;
        } else {
          improving = false;
        }
      }
      
      current_point = next_point;
    }
    
    print(iters);
    
    return current_point;
  };
  
  return corr_fn;
  
}


corrector_fn create_two_point_stepsize_descent(const optim_fn& fn,
                                               const gradient_fn& grad_fn,
                                               double eps_initial_step_size,
                                               double eps_descent_direction,
                                               const double_vector& lower,
                                               const double_vector& upper) {

  corrector_fn corr_fn = [fn, grad_fn, eps_initial_step_size, eps_descent_direction, lower, upper]
                         (evaluated_point starting_point, double_vector ref_point) {
    // Setup
    
    evaluated_point previous_iterate;
    evaluated_point current_iterate;

    double_vector descent_direction = mo_steepest_descent_direction(grad_fn(starting_point));
    double_vector previous_descent_direction = descent_direction;

    double alpha = eps_initial_step_size;
    double scale_factor = 2;
    
    // First: Some line search to find initial alpha
    
    evaluated_point trial_point = starting_point;
    
    do {
      current_iterate = trial_point;
      
      trial_point.dec_space = ensure_boundary(starting_point.dec_space + alpha * normalize(descent_direction),
                                              lower, upper);
      trial_point.obj_space = fn(trial_point.dec_space);
      
      alpha *= scale_factor;
    } while (dominates(trial_point.obj_space, current_iterate.obj_space));
    
    // We increased alpha once too often
    alpha /= scale_factor;
    // Set to starting_point and treat previous step as first iteration
    previous_iterate = starting_point;
    previous_descent_direction = descent_direction;
    
    
    if (current_iterate.dec_space == starting_point.dec_space) {
      // We did not do any progress, so we can stop here
      return starting_point;
    } else {
      descent_direction = mo_steepest_descent_direction(grad_fn(current_iterate));
      
      while (norm(descent_direction) >= eps_descent_direction && abs(alpha) >= eps_initial_step_size) {
        
        double_vector sk = current_iterate.dec_space - previous_iterate.dec_space;
        double_vector yk = descent_direction - previous_descent_direction;
        
        // Stabilized Barzilai-Borwein
        // alpha = min((dot(sk, sk) / dot(sk, -yk)), 2 * alpha);
        double alpha_bb = (dot(sk, sk) / dot(sk, -yk));
        double alpha_pos = norm(sk) / norm(yk);
        alpha = max(alpha_bb, alpha_pos);

        if (alpha == inf) {
          break;
        }
        
        // print_vector(descent_direction);
        // print(alpha);
        // print("");
        
        do {
          trial_point.dec_space = ensure_boundary(current_iterate.dec_space + alpha * descent_direction,
                                        lower, upper);
          trial_point.obj_space = fn(trial_point.dec_space);
          
          alpha /= scale_factor;
        } while (!dominates(trial_point.obj_space, starting_point.obj_space));
        
        // Decreased once too often
        alpha *= scale_factor;
        
        // Update State
        
        previous_iterate = current_iterate;
        previous_descent_direction = descent_direction;
        
        current_iterate = trial_point;
        descent_direction = mo_steepest_descent_direction(grad_fn(current_iterate));
        
      }
    }
    
    return current_iterate;
  };

  return corr_fn;

}

corrector_fn create_adaptive_gradient_descent(const optim_fn& fn,
                                              const gradient_fn& grad_fn,
                                              double eps_initial_step_size,
                                              double eps_descent_direction,
                                              const double_vector& lower,
                                              const double_vector& upper) {
  
  corrector_fn corr_fn = [fn, grad_fn, eps_initial_step_size, eps_descent_direction, lower, upper]
  (evaluated_point starting_point, double_vector ref_point) {
    // Setup
    
    evaluated_point iterate;
    evaluated_point previous_iterate = starting_point;
    
    double_vector descent_direction = mo_steepest_descent_direction(grad_fn(starting_point));
    double_vector previous_descent_direction = descent_direction;

    iterate.dec_space = ensure_boundary(previous_iterate.dec_space + eps_initial_step_size * normalize(descent_direction),
                                                lower, upper);
    iterate.obj_space = fn(iterate.dec_space);
    descent_direction = mo_steepest_descent_direction(grad_fn(iterate));
    
    double theta = inf;
    double lambda = eps_initial_step_size;
    double previous_lambda;
    
    while (norm(descent_direction) > eps_descent_direction) {
      
      lambda = min(
        sqrt(1 + theta) * previous_lambda,
        norm(iterate.dec_space - previous_iterate.dec_space) / (2 * norm(descent_direction - previous_descent_direction))
      );
      
      if (lambda == inf) {
        lambda = 1;
      }

      evaluated_point trial_point;
      
      do {
        trial_point.dec_space = ensure_boundary(iterate.dec_space + lambda * descent_direction, lower, upper);
        trial_point.obj_space = fn(trial_point.dec_space);
        
        lambda /= 2;
      } while (!dominates(trial_point.obj_space, starting_point.obj_space) &&
                lambda * norm(descent_direction) >= eps_initial_step_size);
      
      if (lambda * norm(descent_direction) <= eps_initial_step_size) {
        break;
      }
      
      lambda *= 2;
      
      theta = lambda / previous_lambda;
      
      previous_lambda = lambda;
      previous_iterate = iterate;
      previous_descent_direction = descent_direction;
      
      iterate = trial_point;
      descent_direction = mo_steepest_descent_direction(grad_fn(iterate));
    }
    
    // print(norm(descent_direction));
    
    return iterate;
  };
  
  return corr_fn;
  
}