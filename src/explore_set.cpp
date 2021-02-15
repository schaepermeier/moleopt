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
    // Generic helper variables
    
    evaluated_point previous_point;
    evaluated_point next_point;
    double_vector set_direction;
    
    bool finished = false;
    evaluated_point current_point = starting_point;
    
    double adjusted_explore_set = eps_explore_set;
    
    while (!finished) {
      if (previous_point.dec_space.size() == 0) {
        // estimate direction using single-objective gradient
        // only make a very small step in this case
        set_direction = 0.1 * -normalize(current_gradients[objective]);
      } else {
        set_direction = normalize(current_point.dec_space - previous_point.dec_space);
      }
      
      next_point.dec_space = current_point.dec_space + adjusted_explore_set * set_direction;
      next_point.dec_space = ensure_boundary(next_point.dec_space, lower, upper);
      next_point.obj_space = fn(next_point.dec_space);
      
      while (adjusted_explore_set > eps_explore_set && dominates(current_point.obj_space, next_point.obj_space)) {
        adjusted_explore_set /= 2;
        next_point.dec_space = current_point.dec_space + adjusted_explore_set * set_direction;
        next_point.dec_space = ensure_boundary(next_point.dec_space, lower, upper);
        next_point.obj_space = fn(next_point.dec_space);
      }
      
      double_vector ref_point = next_point.obj_space;
      // ref_point[objective] = current_point.obj_space[objective];
      
      auto descent_point = descent_fn(next_point, ref_point);
      double delta = norm(descent_point.dec_space - next_point.dec_space);
      
      if (!dominates(descent_point.obj_space, current_point.obj_space) &&
          adjusted_explore_set > eps_explore_set &&
          delta >= adjusted_explore_set) {
        adjusted_explore_set /= 2;
        continue;
      } else if (delta <= 0.2 * adjusted_explore_set) {
        adjusted_explore_set *= 1.2;
        adjusted_explore_set = min(adjusted_explore_set, max_explore_set);
      } else {
        adjusted_explore_set /= 1.2;
        adjusted_explore_set = max(adjusted_explore_set, eps_explore_set);
      }
      
      next_point = descent_point;
      
      if (next_point.obj_space[objective] < current_point.obj_space[objective] &&
          !strictly_dominates(next_point.obj_space, current_point.obj_space)) {
          // successfully made step in set
          previous_point = current_point;
        current_point = next_point;
        
        set.insert({current_point.obj_space[0], current_point});
      } else {
        if (strictly_dominates(next_point.obj_space, current_point.obj_space)) {
          // (Presumably) found new attraction basin!
          ridged_points.push_back(next_point);
        } else if (next_point.obj_space[objective] >= current_point.obj_space[objective]) {
          // (Presumably) reached single-objective optimum
          // TODO single-objective descent
        } else {
          // Something else happened that is weird
          print("Something else, WTF");
        }
        // Either way, we are finished here
        finished = true;
      }
    }
  }
  
  // empty vector if finished
  // current point if gone over ridge
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
