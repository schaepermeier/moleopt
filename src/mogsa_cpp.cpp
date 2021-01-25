#include "mogsa_cpp.h"
#include "utils.h"
#include "mo_descent.h"
#include <set>
#include <cmath>

using namespace std;

double_vector lower;
double_vector upper;

double eps_explore_set;
double eps_initial_step_size;
double max_explore_set;

optim_fn fn;
gradient_fn grad_fn;
corrector_fn descend_fn;

std::tuple<evaluated_point, std::vector<evaluated_point>> explore_efficient_set(evaluated_point current_point, int objective) {
  std::vector<double_vector> current_gradients = grad_fn(current_point);

  evaluated_point previous_point;
  evaluated_point next_point;
  double_vector set_direction;
  
  std::vector<evaluated_point> trace;
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
    
    auto descent_point = descend_fn(next_point, ref_point);
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

std::tuple<std::vector<std::map<double, evaluated_point>>,
           std::vector<std::tuple<int, int>>> run_mogsa(optim_fn f, std::vector<double_vector> starting_points, double_vector lower_bounds, double_vector upper_bounds,
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size, double maximum_explore_set) {
  eps_explore_set = epsilon_explore_set;
  eps_initial_step_size = epsilon_initial_step_size;
  max_explore_set = maximum_explore_set;
  
  fn = f;
  lower = lower_bounds;
  upper = upper_bounds;
  
  grad_fn = create_gradient_fn(fn,
                               lower,
                               upper,
                               "twosided",
                               epsilon_gradient);
  
  descend_fn = create_armijo_descent_corrector(fn,
                                               grad_fn,
                                               eps_initial_step_size,
                                               lower,
                                               upper);
  
  std::vector<std::map<double, evaluated_point>> local_sets;
  std::vector<std::tuple<int, int>> set_transitions;
  
  int starting_points_done = 0;
  
  for (double_vector starting_point : starting_points) {
    starting_points_done++;
    print("Starting point No. " + to_string(starting_points_done));
    
    evaluated_point current_point = {
      starting_point,
      fn(starting_point)
    };

    evaluated_point next_point;
    
    std::vector<std::tuple<evaluated_point, int>> points_to_explore;
    
    double_vector ref_point_offset = {0, 0};
    auto descent_point = descend_fn(current_point, current_point.obj_space + ref_point_offset);
    current_point = descent_point;

    points_to_explore.push_back({current_point, -1});

    while(points_to_explore.size() > 0) {
      auto [point_to_explore, origin_set_id] = points_to_explore.back();
      current_point = point_to_explore;
      points_to_explore.pop_back();
      
      // Validate that chosen point does not belong to an already explored set
      
      bool already_explored = false;
      int containing_set = 0;
      
      for (auto& set : local_sets) {
        containing_set++;
        
        double lower_f1 = (*(set.begin())).second.obj_space[0];
        double upper_f2 = (*(set.begin())).second.obj_space[1];
        double upper_f1 = (*(set.rbegin())).second.obj_space[0];
        double lower_f2 = (*(set.rbegin())).second.obj_space[1];
        
        if (current_point.obj_space[0] >= lower_f1 && 
            current_point.obj_space[0] <= upper_f1 && 
            current_point.obj_space[1] >= lower_f2 && 
            current_point.obj_space[1] <= upper_f2) {
          
          std::map<double, evaluated_point>::iterator it_lower = set.lower_bound(current_point.obj_space[0]);
  
          if ((it_lower != set.end()) && (it_lower != set.begin())) {
            auto& right_neighbor = (*it_lower).second;
            auto& left_neighbor = (*(--it_lower)).second;

            if (current_point.obj_space[0] >= left_neighbor.obj_space[0] && 
                current_point.obj_space[0] <= right_neighbor.obj_space[0] && 
                current_point.obj_space[1] >= right_neighbor.obj_space[1] && 
                current_point.obj_space[1] <= left_neighbor.obj_space[1]) {
              
              double norm_to_left = norm(current_point.dec_space - left_neighbor.dec_space);
              double norm_to_right = norm(current_point.dec_space - right_neighbor.dec_space);
              double norm_left_right = norm(left_neighbor.dec_space - right_neighbor.dec_space);

              if ((norm_to_left * norm_to_left + norm_to_right * norm_to_right) <= (norm_left_right * norm_left_right + eps_initial_step_size)) {
              // if (norm_to_left < norm_left_right && norm_to_right < norm_left_right) {
                set.insert(std::pair<double, evaluated_point>(current_point.obj_space[0], current_point));
                already_explored = true;
                break;
              }
            }
          }
        }
      }
      
      if (already_explored) {
        set_transitions.push_back({origin_set_id, containing_set});
        print("Skipping: Set already explored");
        continue;
      }
      
      print("Exploring new set");
      print_vector(current_point.dec_space);
      print("Points left: " + to_string(points_to_explore.size()));
      
      std::map<double, evaluated_point> current_set;
      int current_set_id = local_sets.size() + 1;
      
      set_transitions.push_back({origin_set_id, current_set_id});
      
      for (int obj = 0; obj < 2; obj++) {
        const auto [next_point, trace] = explore_efficient_set(current_point, obj);
        
        if (next_point.dec_space.size() != 0) {
          points_to_explore.push_back({next_point, current_set_id});
        }
        
        for (const auto& eval_point : trace) {
          current_set.insert(std::pair<double, evaluated_point>(eval_point.obj_space[0], eval_point));
        }
      }
      
      local_sets.push_back(current_set);
    }
  }
  
  return {local_sets, set_transitions};
}
