#include "explore_set.h"
#include "vector_utils.h"

tuple<efficient_set, vector<evaluated_point>> explore_efficient_set(
    const evaluated_point& starting_point,
    const optim_fn& fn,
    const gradient_fn& grad_fn,
    const corrector_fn& descent_fn,
    const double_vector& lower,
    const double_vector& upper,
    double eps_explore_set,
    double max_explore_set) {
  // Gradients at starting point, used for
  // first extrapolation in each direction
  vector<double_vector> current_gradients = grad_fn(starting_point);
  
  // Setup efficient set used that is returned later
  efficient_set set;
  set.insert({starting_point.obj_space[0], starting_point});
  
  // Setup vector for points that crossed a ridge
  vector<evaluated_point> ridged_points;
  
  double max_angle_deviation = 45; // maximum angle deviation during set exploration
  double max_step_factor = 2; // max. for difference between two steps in same set
  
  for (int objective = 0; objective < 2; objective++) {

    // Some helper objects
    evaluated_point most_recent;
    evaluated_point second_most_recent;
    evaluated_point predicted;
    evaluated_point corrected;
    
    double_vector set_direction;
    double_vector ref_point;

    double step_size = eps_explore_set;

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
      
      set_direction = normalize(set_direction);
      
      /* (2) Perform prediction along efficient set, check integrity */

      predicted.dec_space = most_recent.dec_space + step_size * set_direction;
      ensure_boundary(predicted.dec_space, lower, upper);
      predicted.obj_space = fn(predicted.dec_space);

      // Check already that new predicted point is better in tracked objective

      if (predicted.obj_space[objective] >= most_recent.obj_space[objective]) {
        // print("Predicted was worse than Most Recent");
        
        if (step_size > eps_explore_set) {
          step_size = max(step_size / max_step_factor, eps_explore_set);
          continue;
        } else {
          if (force_gradient_direction) {
            terminate = true;
          } else {
            force_gradient_direction = true;
          }
          continue;
        }
      }
      
      // TODO terminate early, if dominates(predicted, most_recent)?
      if (dominates(predicted.obj_space, most_recent.obj_space) && step_size > eps_explore_set) {
        print("Reduced early!");
        step_size = max(step_size / max_step_factor, eps_explore_set);
        continue;
      }

      // double_vector ref_point = {inf, inf};
      // ref_point[objective] = most_recent.obj_space[objective];
      ref_point = predicted.obj_space;
      
      double max_descent;
      if (step_size == eps_explore_set && force_gradient_direction) {
        max_descent = inf;
      } else {
        max_descent = step_size;
      }
      
      /* (3) Correct back to locally efficient point, check integrity */
      
      corrected = descent_fn(predicted, ref_point, max_descent);

      double correction_distance = norm(predicted.dec_space - corrected.dec_space);
      
      double angle_to_corrected;
      if (set.size() > 1) {
        angle_to_corrected = angle(most_recent.dec_space - second_most_recent.dec_space,
                                   most_recent.dec_space - corrected.dec_space);
        if (isnan(angle_to_corrected)) {
          // In this case, most likely, the last two steps were almost
          // into the same direction which induced some numerical issue
          angle_to_corrected = 180;
        }
      } else {
        angle_to_corrected = 180;
      }

      if (correction_distance > step_size ||
          angle_to_corrected < (180 - max_angle_deviation)) {

        if (correction_distance > step_size) {
          print("The correction distance was too large " + to_string(correction_distance) +
            "/" + to_string(step_size));
        } else if (angle_to_corrected < (180 - max_angle_deviation)) {
          print("The angle to corrected was too small " + to_string(angle_to_corrected));
        }
        
        if (step_size > eps_explore_set) {
          step_size = max(step_size / 2, eps_explore_set);
          continue;
        } else {
          if (!force_gradient_direction) {
            force_gradient_direction = true;
            continue;
          }
        }
      }
      
      /* (4) If successful until here, determine if we ridged or are in the same set */
      
      if (norm(most_recent.dec_space - corrected.dec_space) > max_explore_set && force_gradient_direction) {
        // Only executed if step_size == eps_explore_set && force_gradient_direction
        ridged_points.push_back(corrected);
        terminate = true;
      } else if (dominates(corrected.obj_space, most_recent.obj_space)) {
        ridged_points.push_back(corrected);
        terminate = true;
      } else {
        set.insert({corrected.obj_space[0], corrected});
        // print(corrected.obj_space - most_recent.obj_space);
        // print(angle_to_corrected);
        // print(correction_distance);
        // print(step_size);
        // print("");
        
        /* (5) Update step size */
        
        if (force_gradient_direction) {
          force_gradient_direction = false;
        } else {
          double angle_deviation = 180 - angle_to_corrected;
          
          double step_size_factor;
          
          if (angle_deviation == 0) {
            // Set to max
            step_size_factor = max_step_factor;
          } else {
            step_size_factor = (max_angle_deviation / 2) / angle_deviation;

            // Keep step size factor reasonable
            step_size_factor = min(step_size_factor, max_step_factor);
            step_size_factor = max(step_size_factor, 1 / max_step_factor);
          }
          
          step_size = step_size * step_size_factor;
          step_size = min(step_size, max_explore_set);
          step_size = max(step_size, eps_explore_set);
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
    double eps_explore_set,
    double max_explore_set) {
  
  explore_set_fn f = [&fn, &grad_fn, &descent_fn, lower, upper, eps_explore_set, max_explore_set]
  (const evaluated_point& starting_point) mutable {
    return explore_efficient_set(
      starting_point,
      fn,
      grad_fn,
      descent_fn,
      lower,
      upper,
      eps_explore_set,
      max_explore_set);
  };
  
  return f;
}
