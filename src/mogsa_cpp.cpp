#include "mogsa_cpp.h"
#include "mo_descent.h"
#include "set_utils.h"
#include "utils.h"
#include <set>
#include <cmath>

double_vector lower;
double_vector upper;

double eps_explore_set;
double eps_initial_step_size;
double max_explore_set;

optim_fn fn;
gradient_fn grad_fn;
corrector_fn descent_fn;

tuple<evaluated_point, vector<evaluated_point>> explore_efficient_set(evaluated_point current_point, int objective) {
  vector<double_vector> current_gradients = grad_fn(current_point);

  evaluated_point previous_point;
  evaluated_point next_point;
  double_vector set_direction;
  
  vector<evaluated_point> trace;
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
      
      trace.push_back(current_point);
    } else {
      if (strictly_dominates(next_point.obj_space, current_point.obj_space)) {
        return {next_point, trace};
      } else if (next_point.obj_space[objective] >= current_point.obj_space[objective]) {
        // print("SO Optimum reached");
        return {{}, trace};
      } else {
        print("Something else lol");
        return {next_point, trace};
      }
      // either hit single-objective optimum, gone over ridge or some error
      finished = true;
    }
  }
  
  // empty vector if finished
  // current point if gone over ridge
  return {{}, trace};
}

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
    optim_fn mo_function,
    gradient_fn gradient_function,
    corrector_fn descent_function,
    vector<double_vector> starting_points,
    double_vector lower_bounds,
    double_vector upper_bounds,
    double epsilon_explore_set,
    double epsilon_initial_step_size,
    double maximum_explore_set) {
  
  /* ========= Setup ========= */
             
  eps_explore_set = epsilon_explore_set;
  max_explore_set = maximum_explore_set;
  
  eps_initial_step_size = epsilon_initial_step_size;
  
  fn = mo_function;
  lower = lower_bounds;
  upper = upper_bounds;
  
  grad_fn = gradient_function;
  descent_fn = descent_function;
  
  set_duplicated_fn already_visited_fn = check_duplicated_set;
  
  /* ========= Mogsa++ Algorithm ========= */
  
  vector<efficient_set> local_sets;
  vector<tuple<int, int>> set_transitions;
  
  int starting_points_done = 0;
  
  for (double_vector starting_point_dec : starting_points) {
    starting_points_done++;
    print("Starting point No. " + to_string(starting_points_done));
    
    evaluated_point starting_point = {
      starting_point_dec,
      fn(starting_point_dec)
    };
    
    double_vector ref_point_offset = {0, 0};
    starting_point = descent_fn(starting_point, starting_point.obj_space + ref_point_offset);

    // The locally efficient points that should be explored
    // during this iteration of the local search.
    // That is the (descended) initial point, and all (descended)
    // points derived from crossed ridges.
    //
    // -1 denotes that a point was a (descended) starting point.
    // Otherwise it denotes the ID of the origin set (before crossing a ridge).
    vector<tuple<evaluated_point, int>> points_to_explore = {
      {starting_point, -1}
    };

    while(points_to_explore.size() > 0) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      points_to_explore.pop_back();
      
      // Validate that chosen point does not belong to an already explored set
      
      int containing_set = already_visited_fn(local_sets, point_to_explore, eps_initial_step_size);
      
      if (containing_set != -1) {
        // This set was already explored before.
        // Log the set transition and continue.
        
        local_sets[containing_set].insert(pair<double, evaluated_point>(point_to_explore.obj_space[0], point_to_explore));
        set_transitions.push_back({origin_set_id, containing_set});
        print("Skipping: Set already explored");
      } else {
        // This point belongs to an unexplored set
        
        print("Exploring new set");
        print_vector(point_to_explore.dec_space);
        print("Points left: " + to_string(points_to_explore.size()));
        
        // Create new local set
        efficient_set current_set;
        int current_set_id = local_sets.size() + 1;
        
        // Log set transition
        set_transitions.push_back({origin_set_id, current_set_id});
        
        // TODO make explore_efficient_set responsible for tracking objectives
        
        for (int obj = 0; obj < 2; obj++) {
          const auto [next_point, trace] = explore_efficient_set(point_to_explore, obj);
          
          if (next_point.dec_space.size() != 0) {
            // We crossed a ridge and got a (descended) point
            // that should be explored in this iteration
            points_to_explore.push_back({next_point, current_set_id});
          }
          
          // Add all newly discovered points to the current set
          for (const auto& eval_point : trace) {
            current_set.insert(pair<double, evaluated_point>(eval_point.obj_space[0], eval_point));
          }
        }
        
        local_sets.push_back(current_set);
      }
      
    }
  }
  
  return {local_sets, set_transitions};
}
