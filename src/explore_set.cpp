#include "explore_set.h"
#include "vector_utils.h"

tuple<efficient_set, vector<evaluated_point>> explore_efficient_set(
    evaluated_point starting_point,
    optim_fn fn,
    gradient_fn grad_fn,
    corrector_fn descent_fn,
    double_vector lower,
    double_vector upper,
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

    double step_size = eps_explore_set;

    bool terminate = false;
    bool force_gradient_direction = false;

    while (!terminate) {
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

      double_vector set_direction;

      if (set.size() == 1 || force_gradient_direction) {
        set_direction = project_feasible_direction(-current_gradients[objective], starting_point.dec_space, lower, upper);
        set_direction = normalize(set_direction);
      } else {
        set_direction = project_feasible_direction(most_recent.dec_space - second_most_recent.dec_space,
                                                   most_recent.dec_space, lower, upper);
        set_direction = normalize(set_direction);
      }
      
      if (norm(set_direction) == 0) {
        break;
      }

      predicted.dec_space = ensure_boundary(most_recent.dec_space + step_size * set_direction, lower, upper);
      predicted.obj_space = fn(predicted.dec_space);

      // Check already that new predicted point is better in tracked objective

      if (predicted.obj_space[objective] >= most_recent.obj_space[objective]) {
        print("Predicted was worse than Most Recent");
        
        if (step_size > eps_explore_set) {
          step_size = max(step_size / 2, eps_explore_set);
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

      // double_vector ref_point = {inf, inf};
      // ref_point[objective] = most_recent.obj_space[objective];
      double_vector ref_point = predicted.obj_space;
      
      double max_descent;
      if (step_size == eps_explore_set && force_gradient_direction) {
        max_descent = inf;
      } else {
        max_descent = step_size;
      }
      
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
      
      if (norm(most_recent.dec_space - corrected.dec_space) > max_explore_set && set.size() != 1 && force_gradient_direction) {
        // Only executed if step_size == eps_explore_set && force_gradient_direction
        ridged_points.push_back(corrected);
        terminate = true;
      } else
        if (strictly_dominates(corrected.obj_space, most_recent.obj_space)) {
        ridged_points.push_back(corrected);
        terminate = true;
      } else if (dominates(corrected.obj_space, most_recent.obj_space)) {
        // We hit a weakly dominated area!
        print("Bruh");
        terminate = true;
      } else {
        set.insert({corrected.obj_space[0], corrected});
        // print(corrected.obj_space - most_recent.obj_space);
        // print(angle_to_corrected);
        // print(correction_distance);
        // print(step_size);
        // print("");
        
        if (force_gradient_direction) {
          force_gradient_direction = false;
        } else {
          double angle_deviation = 180 - angle_to_corrected;
          
          double step_size_factor;
          
          if (angle_deviation == 0) {
            // Set to max
            step_size_factor = 2;
          } else {
            step_size_factor = max_angle_deviation / angle_deviation;
            
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
  
  // Explored local set and all "ridged" points
  return {set, ridged_points};
}

explore_set_fn get_explore_set_fn(
    optim_fn fn,
    gradient_fn grad_fn,
    corrector_fn descent_fn,
    double_vector lower,
    double_vector upper,
    double eps_explore_set,
    double max_explore_set) {
  
  explore_set_fn f = [fn, grad_fn, descent_fn, lower, upper, eps_explore_set, max_explore_set](evaluated_point starting_point) {
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
