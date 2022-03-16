#include "mole.h"
#include "mo_descent.h"
#include "set_utils.h"
#include "explore_set.h"
#include "utils.h"
#include <set>
#include <cmath>

bool evaluate_starting_points = false;

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mole(
    const optim_fn& mo_function,
    const gradient_fn& gradient_function,
    const corrector_fn& descent_function,
    const explore_set_fn& explore_set_function,
    const vector<double_vector>& starting_points,
    const double_vector& lower_bounds,
    const double_vector& upper_bounds,
    int max_local_sets,
    double explore_step_min,
    int refine_after_nstarts,
    double refine_hv_target,
    long max_budget,
    long* used_budget) {
  
  /* ========= Setup ========= */
  
  set_duplicated_fn already_visited_fn = check_duplicated_set;
  
  /* ========= MOLE Algorithm ========= */
  
  vector<efficient_set> local_sets;
  vector<tuple<int, int>> set_transitions;
  
  set<double_vector> nondominated_points;
  
  bool budget_depleted = false;
  
  vector<evaluated_point> starting_points_evaluated;
  
  for (double_vector starting_point_dec : starting_points) {
    if (evaluate_starting_points) {
      starting_points_evaluated.push_back({
        starting_point_dec,
        mo_function(starting_point_dec)
      });
    } else {
      starting_points_evaluated.push_back({
        starting_point_dec,
        {inf, inf}
      });
    }
  }
  
  starting_points_evaluated = sort_by_nondominated(starting_points_evaluated);
  bool did_hv_refinement = false;
  
  for (int point_index = 0; point_index < starting_points_evaluated.size(); point_index++) {
    evaluated_point starting_point = starting_points_evaluated[point_index];
    print_info("Starting point No. " + to_string(point_index + 1));
    // print_info("Evals: ");
    // print_info(*(used_budget));
    // print_info(starting_point.dec_space);
    // print_info(starting_point.obj_space);
    
    bool new_nondominated_set = false;
    
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
    
    while (points_to_explore.size() > 0 && local_sets.size() < max_local_sets) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      points_to_explore.pop_back();
      
      if (!evaluate_starting_points && point_to_explore.obj_space[0] == inf) {
        point_to_explore.obj_space = mo_function(point_to_explore.dec_space);
      }
      
      if (point_to_explore.obj_space[0] == inf) {
        // Terminate if budget is empty, i.e.
        // fn == inf
        budget_depleted = true;
        break;
      }
      
      // Validate that chosen point does not belong to an already explored set
      
      int containing_set = already_visited_fn(local_sets, point_to_explore, explore_step_min);
      
      if (containing_set != -1) {
        // This set was already explored before.
        // Log the set transition and continue.
        
        print_info("Skipping: Set already explored");
        
        // local_sets[containing_set].insert(pair<double, evaluated_point>(point_to_explore.obj_space[0], point_to_explore));
        insert_into_set(local_sets[containing_set], point_to_explore);
        insert_nondominated(nondominated_points, point_to_explore.obj_space);
        
        // Avoid loops for individual sets
        if (origin_set_id != containing_set) {
          set_transitions.push_back({origin_set_id, containing_set});
        }
      } else {
        // This point belongs to an unexplored set
        
        print_info("Exploring new set (No. " + to_string(local_sets.size() + 1) + ")");
        print(point_to_explore.dec_space);
        print("Points left: " + to_string(points_to_explore.size()));
        
        if (point_to_explore.obj_space[0] == inf) {
          // Terminate if budget is empty, i.e.
          // fn == inf
          print_info("Terminating: Budget used up!");
          break;
        }
        
        // Explore the new local set
        auto [set, ridged_points] = explore_set_function(point_to_explore);
        
        // Store the new local set
        int set_id = local_sets.size();
        local_sets.push_back(set);
        
        for (auto& [f1_val, point] : set) {
          bool insert_successful = insert_nondominated(nondominated_points, point.obj_space);
          new_nondominated_set = new_nondominated_set && insert_successful;
        }

        // Log set transition
        set_transitions.push_back({origin_set_id, set_id});
        
        // Queue new points to explore (if any)
        for (evaluated_point point : ridged_points) {
          points_to_explore.push_back({point, set_id});
        }
      }
    }
    
    /* ========= Post-Processing for HV Maximization ========= */
    
    if (refine_hv_target > 0) {
      if (evaluate_starting_points && (!did_hv_refinement || new_nondominated_set)) {
        
        bool all_dominated = true;
        for (int upcoming_index = point_index + 1; upcoming_index < starting_points_evaluated.size(); upcoming_index++) {
          all_dominated = all_dominated &&
            point_dominated_by_set(starting_points_evaluated[upcoming_index].obj_space, nondominated_points);
        }
        
        if (all_dominated) {
          print_info("refining after " + to_string(point_index + 1) + " points");
          refine_sets(local_sets, refine_hv_target, mo_function, descent_function, nondominated_points);
          
          did_hv_refinement = true;
        }
      } else if ((point_index + 1 >= refine_after_nstarts) && (!did_hv_refinement || new_nondominated_set)) {
        refine_sets(local_sets, refine_hv_target, mo_function, descent_function, nondominated_points);
        
        did_hv_refinement = true;
      }
    }
    
    if (budget_depleted) break;
  }
  
  print("Size of nondominated set: " + to_string(nondominated_points.size()));
  
  return {local_sets, set_transitions};
}

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mole(
    const optim_fn& mo_function,
    const vector<double_vector>& starting_points,
    const double_vector& lower,
    const double_vector& upper,
    int max_local_sets = 1000,
    double epsilon_gradient = 1e-8,
    double descent_direction_min = 1e-8,
    double descent_step_min = 1e-6,
    double descent_step_max = 1e-1,
    double descent_scale_factor = 2,
    double descent_armijo_factor = 1e-4,
    int descent_history_size = 100,
    int descent_max_iter = 1000,
    double explore_step_min = 1e-4,
    double explore_step_max = 1e-1,
    double explore_angle_max = 45,
    double explore_scale_factor = 2,
    int refine_after_nstarts = 10,
    double refine_hv_target = 2e-5) {
  
  /* ========= Setup and run MOLE ========= */
  
  // Create Gradient of fn
  
  gradient_fn gradient_function = create_gradient_fn(mo_function,
                                                     lower,
                                                     upper,
                                                     "twosided",
                                                     epsilon_gradient);
  
  // Create the descent function
  
  corrector_fn descent_function;
  
  // descent_function = create_two_point_stepsize_descent(mo_function,
  //                                                      gradient_function,
  //                                                      lower,
  //                                                      upper,
  //                                                      descent_direction_min,
  //                                                      descent_step_min,
  //                                                      descent_step_max,
  //                                                      descent_scale_factor,
  //                                                      descent_armijo_factor,
  //                                                      descent_history_size,
  //                                                      descent_max_iter);

  descent_function = create_slow_mo_descent(mo_function,
                                           gradient_function,
                                           lower,
                                           upper,
                                           descent_direction_min,
                                           descent_step_min,
                                           descent_step_max,
                                           descent_scale_factor,
                                           descent_armijo_factor,
                                           descent_history_size,
                                           descent_max_iter);

  // Create the explore_set function
  
  explore_set_fn explore_set_function = get_explore_set_fn(
    mo_function,
    gradient_function,
    descent_function,
    lower,
    upper,
    explore_step_min,
    explore_step_max,
    explore_angle_max,
    explore_scale_factor);
  
  // Run MOLE
  
  long max_budget = -1;
  long used_budget = 0;
  
  return run_mole(
    mo_function,
    gradient_function,
    descent_function,
    explore_set_function,
    starting_points,
    lower,
    upper,
    max_local_sets,
    explore_step_min,
    refine_after_nstarts,
    refine_hv_target,
    max_budget,
    (&used_budget));
  
}
