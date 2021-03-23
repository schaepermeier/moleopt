#include "mogsa_cpp.h"
#include "mo_descent.h"
#include "set_utils.h"
#include "explore_set.h"
#include "utils.h"
#include <set>
#include <cmath>

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
    const optim_fn& mo_function,
    const gradient_fn& gradient_function,
    const corrector_fn& descent_function,
    const explore_set_fn& explore_set_function,
    const vector<double_vector>& starting_points,
    const double_vector& lower_bounds,
    const double_vector& upper_bounds,
    double epsilon_explore_set,
    double epsilon_initial_step_size,
    double maximum_explore_set) {
  
  /* ========= Setup ========= */
  
  set_duplicated_fn already_visited_fn = check_duplicated_set;
  
  /* ========= Mogsa++ Algorithm ========= */
  
  vector<efficient_set> local_sets;
  vector<tuple<int, int>> set_transitions;
  
  set<double_vector> nondominated_points;

  int starting_points_done = 0;
  
  bool budget_depleted = false;
  
  for (double_vector starting_point_dec : starting_points) {
    starting_points_done++;
    bool nondom_set_in_iter = false;
    print_info("Starting point No. " + to_string(starting_points_done));
    
    evaluated_point starting_point = {
      starting_point_dec,
      mo_function(starting_point_dec)
    };
    
    if (budget_depleted) {
      // Terminate if budget is empty, i.e.
      // fn == inf
      print_info("Terminating: Budget used up!");
      break;
    }
    
    double_vector ref_point_offset = {0, 0};
    // double_vector ref_point_offset = {inf, inf};
    starting_point = descent_function(starting_point, starting_point.obj_space + ref_point_offset, inf);

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

    while (points_to_explore.size() > 0 && local_sets.size() < 1000) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      points_to_explore.pop_back();
      
      if (point_to_explore.obj_space[0] == inf) {
        // Terminate if budget is empty, i.e.
        // fn == inf
        budget_depleted = true;
        break;
      }
      
      // Validate that chosen point does not belong to an already explored set
      
      int containing_set = already_visited_fn(local_sets, point_to_explore, epsilon_explore_set);
      
      if (containing_set != -1) {
        // This set was already explored before.
        // Log the set transition and continue.
        
        print_info("Skipping: Set already explored");
        
        // local_sets[containing_set].insert(pair<double, evaluated_point>(point_to_explore.obj_space[0], point_to_explore));
        insert_into_set(local_sets[containing_set], point_to_explore);
        insert_nondominated(nondominated_points, point_to_explore.obj_space);
        
        set_transitions.push_back({origin_set_id, containing_set});
      } else {
        // This point belongs to an unexplored set
        
        print_info("Exploring new set (No. " + to_string(local_sets.size() + 1) + ")");
        print(point_to_explore.dec_space);
        print("Points left: " + to_string(points_to_explore.size()));
        
        // Explore the new local set
        auto [set, ridged_points] = explore_set_function(point_to_explore);
        
        // Store the new local set
        int set_id = local_sets.size();
        local_sets.push_back(set);
        
        int nondom_before = nondominated_points.size();
        
        print_info("Inserting points nondominated");
        
        for (auto& [f1_val, point] : set) {
          insert_nondominated(nondominated_points, point.obj_space);
        }
        
        if (nondominated_points.size() > nondom_before) {
         nondom_set_in_iter = true;
        }
        
        // Log set transition
        set_transitions.push_back({origin_set_id, set_id});
        
        // Queue new points to explore (if any)
        for (evaluated_point point : ridged_points) {
          points_to_explore.push_back({point, set_id});
        }
        
        print_info("Finished exploration");
      }
      
      if (local_sets.size() >= 1000) {
        break;
      }
    }
    
    if (starting_points_done == 10 || (starting_points_done > 10 && nondom_set_in_iter)) {
      refine_sets(local_sets, 2e-5, mo_function, descent_function, nondominated_points);
      break;
    }
  }
  
  if (!budget_depleted) {
    refine_sets(local_sets, 1e-5, mo_function, descent_function, nondominated_points);
  }
  
  print(nondominated_points.size());
  
  return {local_sets, set_transitions};
}

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
    const optim_fn& mo_function,
    const vector<double_vector>& starting_points,
    const double_vector& lower,
    const double_vector& upper,
    double epsilon_gradient,
    double epsilon_explore_set,
    double epsilon_initial_step_size,
    double max_explore_set) {
  
  /* ========= Setup and run Mogsa ========= */
  
  // Create Gradient of fn
  
  gradient_fn gradient_function = create_gradient_fn(mo_function,
                                                     lower,
                                                     upper,
                                                     "twosided",
                                                     epsilon_gradient);
  
  // Create the descent function
  
  corrector_fn descent_function;
  
  descent_function = create_two_point_stepsize_descent(mo_function,
                                                       gradient_function,
                                                       epsilon_initial_step_size,
                                                       1e-8,
                                                       lower,
                                                       upper);
  
  // Create the explore_set function
  
  explore_set_fn explore_set_function = get_explore_set_fn(
    mo_function,
    gradient_function,
    descent_function,
    lower,
    upper,
    epsilon_explore_set,
    max_explore_set);
  
  // Run Mogsa
  
  return run_mogsa(
    mo_function,
    gradient_function,
    descent_function,
    explore_set_function,
    starting_points,
    lower,
    upper,
    epsilon_explore_set,
    epsilon_initial_step_size,
    max_explore_set);
  
}

