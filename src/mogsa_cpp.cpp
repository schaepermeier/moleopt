#include "mogsa_cpp.h"
#include "mo_descent.h"
#include "set_utils.h"
#include "utils.h"
#include <set>
#include <cmath>

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
    optim_fn mo_function,
    gradient_fn gradient_function,
    corrector_fn descent_function,
    explore_set_fn explore_set_function,
    vector<double_vector> starting_points,
    double_vector lower_bounds,
    double_vector upper_bounds,
    double epsilon_explore_set,
    double epsilon_initial_step_size,
    double maximum_explore_set) {
  
  /* ========= Setup ========= */
  
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
      mo_function(starting_point_dec)
    };
    
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

    while (points_to_explore.size() > 0) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      points_to_explore.pop_back();
      
      // Validate that chosen point does not belong to an already explored set
      
      int containing_set = already_visited_fn(local_sets, point_to_explore, epsilon_explore_set);
      
      if (containing_set != -1) {
        // This set was already explored before.
        // Log the set transition and continue.
        
        print("Skipping: Set already explored");
        
        local_sets[containing_set].insert(pair<double, evaluated_point>(point_to_explore.obj_space[0], point_to_explore));
        set_transitions.push_back({origin_set_id, containing_set});
      } else {
        // This point belongs to an unexplored set
        
        print("Exploring new set");
        print_vector(point_to_explore.dec_space);
        print("Points left: " + to_string(points_to_explore.size()));
        
        // Explore the new local set
        auto [set, ridged_points] = explore_set_function(point_to_explore);
        
        // Store the new local set
        int set_id = local_sets.size();
        local_sets.push_back(set);
        
        // Log set transition
        set_transitions.push_back({origin_set_id, set_id});
        
        // Queue new points to explore (if any)
        for (evaluated_point point : ridged_points) {
          points_to_explore.push_back({point, set_id});
        }
      }
    }
  }
  
  return {local_sets, set_transitions};
}
