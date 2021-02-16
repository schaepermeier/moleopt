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
  
  for (int objective = 0; objective < 2; objective++) {

    // Some helper objects
    evaluated_point most_recent;
    evaluated_point second_most_recent;
    evaluated_point predicted;
    evaluated_point corrected;

    double step_size = eps_explore_set;

    bool terminate = false;

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

      // print(step_size);
      
      double_vector set_direction;

      if (set.size() == 1) {
        set_direction = -normalize(current_gradients[objective]);
      } else {
        set_direction = normalize(most_recent.dec_space - second_most_recent.dec_space);
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
          break;
        }
      }

      double_vector ref_point = {inf, inf};
      ref_point[objective] = most_recent.obj_space[objective];
      corrected = descent_fn(predicted, ref_point);

      double correction_distance = norm(predicted.dec_space - corrected.dec_space);
      
      double angle_to_corrected;
      if (set.size() > 1) {
        angle_to_corrected = angle(most_recent.dec_space - second_most_recent.dec_space,
                                   most_recent.dec_space - corrected.dec_space);
      } else {
        angle_to_corrected = 180;
      }

      if (correction_distance > step_size ||
          angle_to_corrected < 170) {
        print("The angle to corrected was too small " + to_string(angle_to_corrected) +
              " or correction distance too large " + to_string(correction_distance));
        
        if (step_size > eps_explore_set) {
          step_size = max(step_size / 2, eps_explore_set);
          continue;
        }
      }

      if (dominates(corrected.obj_space, most_recent.obj_space)) {
        ridged_points.push_back(corrected);
        terminate = true;
      } else {
        set.insert({corrected.obj_space[0], corrected});
        step_size = min(step_size * 2, max_explore_set);
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
