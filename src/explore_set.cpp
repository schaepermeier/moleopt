#include "explore_set.h"
#include "vector_utils.h"

tuple<efficient_set, vector<evaluated_point>> explore_efficient_set(
    const evaluated_point& starting_point,
    const optim_fn& fn,
    const gradient_fn& grad_fn,
    const corrector_fn& descent_fn,
    const double_vector& lower,
    const double_vector& upper,
    double explore_step_min,
    double explore_step_max,
    double explore_angle_max,
    double explore_scale_factor) {
  // Gradients at starting point, used for
  // first extrapolation in each direction
  vector<double_vector> current_gradients = grad_fn(starting_point);
  
  // Setup efficient set used that is returned later
  efficient_set set;
  set.insert({starting_point.obj_space[0], starting_point});
  
  // Setup vector for points that crossed a ridge
  vector<evaluated_point> ridged_points;

  for (int objective = 0; objective < 2; objective++) {

    // Some helper objects
    evaluated_point most_recent;
    evaluated_point second_most_recent;
    evaluated_point predicted;
    evaluated_point corrected;
    
    double_vector set_direction;
    double_vector ref_point;

    double step_size = explore_step_min;

    bool terminate = false;
    bool force_gradient_direction = false;

    while (!terminate) {
      /* (1) Determine search direction along efficient set */
      
      // If objective == 0, the most recent points are at the beginning of set
      // If objective == 1, the most recent points are at the end of set
      
      if (objective == 0) {
        most_recent = (*(set.begin())).second;

        if (set.size() > 1) {
          second_most_recent = (*(++set.begin())).second;
        }
      } else {
        most_recent = (*(set.rbegin())).second;

        if (set.size() > 1) {
          second_most_recent = (*(++set.rbegin())).second;
        }
      }
      
      if (set.size() == 1) {
        force_gradient_direction = true;
      }

      if (force_gradient_direction) {
        set_direction = -current_gradients[objective];
      } else {
        set_direction = most_recent.dec_space - second_most_recent.dec_space;
      }
      
      project_feasible_direction(set_direction, most_recent.dec_space, lower, upper);
      
      if (norm(set_direction) == 0) {
        break;
      }

      /* (2) Perform prediction along efficient set, check integrity */

      predicted.dec_space = most_recent.dec_space + step_size * set_direction / norm(set_direction);
      ensure_boundary(predicted.dec_space, lower, upper);
      predicted.obj_space = fn(predicted.dec_space);

      // Check already that new predicted point is better in tracked objective
      if (predicted.obj_space[objective] >= most_recent.obj_space[objective]) {
        print("Predicted was worse than Most Recent");
        
        if (step_size > explore_step_min) {
          step_size = max(step_size / explore_scale_factor, explore_step_min);
        } else {
          if (force_gradient_direction) {
            terminate = true;
          } else {
            force_gradient_direction = true;
          }
        }
        
        continue;
      }
      
      // Reduce step size without descent, if we already dominate
      if (dominates(predicted.obj_space, most_recent.obj_space) && step_size > explore_step_min) {
        print("Reduced early!");
        step_size = max(step_size / explore_scale_factor, explore_step_min);
        continue;
      }

      // double_vector ref_point = {inf, inf};
      // ref_point[objective] = most_recent.obj_space[objective];
      ref_point = predicted.obj_space;
      
      double max_descent;
      if (step_size == explore_step_min && force_gradient_direction) {
        max_descent = inf;
      } else {
        max_descent = step_size;
      }
      
      /* (3) Correct back to locally efficient point, check integrity */
      
      corrected = descent_fn(predicted, ref_point, max_descent);

      double correction_distance = norm(predicted.dec_space - corrected.dec_space);
      
      double angle_to_corrected;
      if (set.size() > 1) {
        angle_to_corrected = angle(second_most_recent.dec_space - most_recent.dec_space,
                                   most_recent.dec_space - corrected.dec_space);
        if (isnan(angle_to_corrected)) {
          // In this case, most likely, the last two steps were almost
          // into the same direction which induced some numerical issue
          angle_to_corrected = 0;
        }
      } else {
        angle_to_corrected = 0;
      }

      if (correction_distance > step_size ||
          angle_to_corrected > explore_angle_max) {

        if (correction_distance > step_size) {
          print("The correction distance was too large " + to_string(correction_distance) +
            "/" + to_string(step_size));
        } else if (angle_to_corrected > explore_angle_max) {
          print("The angle deviation was too big: " + to_string(angle_to_corrected));
        }
        
        if (step_size > explore_step_min) {
          step_size = max(step_size / explore_scale_factor, explore_step_min);
          continue;
        } else {
          if (!force_gradient_direction) {
            force_gradient_direction = true;
            continue;
          }
        }
      }
      
      /* (4) If successful until here, determine if we ridged or are in the same set */
      
      if (dominates(corrected.obj_space, most_recent.obj_space) ||
         (norm(most_recent.dec_space - corrected.dec_space) > explore_step_max && force_gradient_direction)) {
        // Only executed if step_size == explore_step_min && force_gradient_direction
        ridged_points.push_back(corrected);
        terminate = true;
      } else {
        set.insert({corrected.obj_space[0], corrected});
        
        /* (5) Update step size */
        
        if (force_gradient_direction) {
          force_gradient_direction = false;
        } else {
          double step_size_factor;
          
          if (angle_to_corrected == 0) {
            // Set to max
            step_size_factor = explore_scale_factor;
          } else {
            step_size_factor = (explore_angle_max / 2) / angle_to_corrected;

            // Keep step size factor reasonable
            step_size_factor = min(step_size_factor, explore_scale_factor);
            step_size_factor = max(step_size_factor, 1 / explore_scale_factor);
          }
          
          step_size = step_size * step_size_factor;
          step_size = min(step_size, explore_step_max);
          step_size = max(step_size, explore_step_min);
        }
      }
    }
  }
  
  print("Final set size: " + to_string(set.size()) + ", ridged points: " + to_string(ridged_points.size()));
  
  // Explored local set and all "ridged" points
  return {set, ridged_points};
}

explore_set_fn get_explore_set_fn(
    const optim_fn& fn,
    const gradient_fn& grad_fn,
    const corrector_fn& descent_fn,
    const double_vector& lower,
    const double_vector& upper,
    double explore_step_min,
    double explore_step_max,
    double explore_angle_max,
    double explore_scale_factor) {
  
  explore_set_fn f = [&fn, &grad_fn, &descent_fn, lower, upper,
                      explore_step_min, explore_step_max,
                      explore_angle_max, explore_scale_factor]
  (const evaluated_point& starting_point) mutable {
    return explore_efficient_set(
      starting_point,
      fn,
      grad_fn,
      descent_fn,
      lower,
      upper,
      explore_step_min,
      explore_step_max,
      explore_angle_max,
      explore_scale_factor);
  };
  
  return f;
}
