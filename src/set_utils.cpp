#include "set_utils.h"
#include <set>

bool inbounds(const double_vector& p,
              const double_vector& lower,
              const double_vector& upper) {
  for (int i = 0; i < p.size(); i++) {
    if ((p[i] < lower[i]) || (p[i] > upper[i])) {
      return false;
    }
  }
  
  return true;
}

int check_duplicated_set(const vector<efficient_set>& local_sets,
                         const evaluated_point& new_point,
                         double epsilon) {
  int containing_set = -1;
  
  for (auto& set : local_sets) {
    containing_set++;
    
    // Check whether the new point is in the "box" of the
    // objective values covered by the given set
    
    double lower_f1 = (*(set.begin())).second.obj_space[0];
    double upper_f2 = (*(set.begin())).second.obj_space[1];
    double upper_f1 = (*(set.rbegin())).second.obj_space[0];
    double lower_f2 = (*(set.rbegin())).second.obj_space[1];
    
    if (inbounds(new_point.obj_space, {lower_f1, lower_f2}, {upper_f1, upper_f2})) {
      
      // Find the two points of the efficient set
      // that the "new point" has to be between, if it is in this set
      
      map<double, evaluated_point>::const_iterator it_lower = set.lower_bound(new_point.obj_space[0]);
      
      if ((it_lower != set.end()) && (it_lower != set.begin())) {
        auto& right_neighbor = (*it_lower).second;
        auto& left_neighbor = (*(--it_lower)).second;
        
        // If our point really is "between" the objective values of two points of that set
        
        // if (new_point.obj_space[0] >= left_neighbor.obj_space[0] && 
        //     new_point.obj_space[0] <= right_neighbor.obj_space[0] && 
        //     new_point.obj_space[1] >= right_neighbor.obj_space[1] && 
        //     new_point.obj_space[1] <= left_neighbor.obj_space[1]) {
        if (inbounds(new_point.obj_space, {left_neighbor.obj_space[0], right_neighbor.obj_space[1]},
                                          {right_neighbor.obj_space[0], left_neighbor.obj_space[1]})) {
          
          double norm_to_left = norm(new_point.dec_space - left_neighbor.dec_space);
          double norm_to_right = norm(new_point.dec_space - right_neighbor.dec_space);
          double norm_left_right = norm(left_neighbor.dec_space - right_neighbor.dec_space);
          
          // Check that they have at most an angle of 90 degrees in decision space
          
          if ((norm_to_left * norm_to_left + norm_to_right * norm_to_right) <= (norm_left_right * norm_left_right + epsilon)) {
            return containing_set;
          }
        }
      }
    }
    
    // Additionally, check whether the selected point is epsilon close to
    // any logged point in any set.
    for (auto& [f1_val, point] : set) {
      if (norm(new_point.dec_space - point.dec_space) < epsilon) {
        return containing_set;
      }
    }
  }
  
  return -1;
}

void insert_into_set(efficient_set& set,
                     const evaluated_point& new_point) {
  double_vector to_delete;
  
  for (auto& [f1_val, point] : set) {
    if (dominates(new_point.obj_space, point.obj_space)) {
      to_delete.push_back(f1_val);
    } else if (dominates(point.obj_space, new_point.obj_space)) {
      // Don't change set, if new point is dominated by it
      return;
    }
  }
  
  for (double v : to_delete) {
    set.erase(v);
  }
  
  set.insert({new_point.obj_space[0], new_point});
}

bool sort_by_first(const tuple<double, int, evaluated_point, evaluated_point>& a,  
               const tuple<double, int, evaluated_point, evaluated_point>& b) { 
  return (get<0>(a) < get<0>(b)); 
}

void refine_sets(vector<efficient_set>& sets,
                 double rel_hv_target,
                 optim_fn fn,
                 corrector_fn descent_fn) {
  efficient_set nondominated_points;
  
  for (auto& set : sets) {
    for (auto& [f1_val, point] : set) {
      insert_into_set(nondominated_points, point);
    }
  }
  
  // Get objective space coordinates of (approximate) nadir and ideal points
  
  evaluated_point opt_f1 = (*(nondominated_points.begin())).second;
  evaluated_point opt_f2 = (*(nondominated_points.rbegin())).second;

  double lower_f1 = opt_f1.obj_space[0];
  double upper_f2 = opt_f1.obj_space[1];
  double upper_f1 = opt_f2.obj_space[0];
  double lower_f2 = opt_f2.obj_space[1];

  // Maximal HV in this approximated objective space (for normalization)

  double max_hv = (upper_f1 - lower_f1) * (upper_f2 - lower_f2);

  double total_hv_potential = 0.0;

  // Iterate over each set.
  // If one of two neighboring points is nondominated, compute potential HV gain between them

  std::vector<std::tuple<double, int, evaluated_point, evaluated_point>> potential_pairs;

  for (int set_id = 0; set_id < sets.size(); set_id++) {
    auto it = sets[set_id].begin();
    
    evaluated_point prev; // "Previous" point along set
    evaluated_point current; // "Current" point along set
    
    bool prev_nondom; // Does the previous point along the set belong to the non-dominated ones?
    bool current_nondom; // Does the current point along the set belong to the non-dominated ones?
    
    current = it->second;

    auto nondom_it = (nondominated_points.find(current.obj_space[0]));
    current_nondom = (nondom_it != nondominated_points.end()) && ((*nondom_it).second.dec_space == current.dec_space);
    
    it++;

    for (; it != sets[set_id].end(); it++) {
      prev = current;
      prev_nondom = current_nondom;
      
      current = it->second;

      nondom_it = (nondominated_points.find(current.obj_space[0]));
      current_nondom = (nondom_it != nondominated_points.end()) && ((*nondom_it).second.dec_space == current.dec_space);

      if (prev_nondom || current_nondom) {
        double_vector improvement = current.obj_space - prev.obj_space;
        
        double hv_potential = (current.obj_space[0] - prev.obj_space[0]) *
                              (prev.obj_space[1] - current.obj_space[1]);
        
        total_hv_potential += hv_potential;
        
        evaluated_point left(prev);
        evaluated_point right(current);
        
        potential_pairs.push_back({hv_potential, set_id, left, right});
      }
    }
  }
  
  while (total_hv_potential / max_hv > rel_hv_target && potential_pairs.size() > 0) {
    sort(potential_pairs.begin(), potential_pairs.end(), sort_by_first);
    
    auto [hv_potential, set_id, left, right] = potential_pairs.back();

    total_hv_potential -= hv_potential;
    potential_pairs.pop_back();
    
    evaluated_point new_point;
    new_point.dec_space = (left.dec_space + right.dec_space) / 2;
    new_point.obj_space = fn(new_point.dec_space);

    if (inbounds(new_point.obj_space, {left.obj_space[0], right.obj_space[1]},
                                      {right.obj_space[0], left.obj_space[1]})) {
      new_point = descent_fn(new_point, new_point.obj_space, inf);

      insert_into_set(sets[set_id], new_point);
      
      double_vector improvement;
      double hv_potential;
      
      // Add {left, new_point}
      
      hv_potential = (new_point.obj_space[0] - left.obj_space[0]) *
        (left.obj_space[1] - new_point.obj_space[1]);
      
      total_hv_potential += hv_potential;
      
      potential_pairs.push_back({hv_potential, set_id, left, new_point});
      
      // Add {new_point, right}
      
      hv_potential = (right.obj_space[0] - new_point.obj_space[0]) *
        (new_point.obj_space[1] - right.obj_space[1]);
      
      total_hv_potential += hv_potential;
      
      potential_pairs.push_back({hv_potential, set_id, new_point, right});
      
      print(total_hv_potential / max_hv);
    } else {
      // print("Sad HV Noises");
    }
  
  }
  
  print(total_hv_potential / max_hv);
}
