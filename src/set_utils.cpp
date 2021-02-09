#include "set_utils.h"

int check_duplicated_set(std::vector<efficient_set> local_sets,
                         evaluated_point new_point,
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
    
    if (new_point.obj_space[0] >= lower_f1 && 
        new_point.obj_space[0] <= upper_f1 && 
        new_point.obj_space[1] >= lower_f2 && 
        new_point.obj_space[1] <= upper_f2) {
      
      // Find the two points of the efficient set
      // that the "new point" has to be between, if it is in this set
      
      std::map<double, evaluated_point>::iterator it_lower = set.lower_bound(new_point.obj_space[0]);
      
      if ((it_lower != set.end()) && (it_lower != set.begin())) {
        auto& right_neighbor = (*it_lower).second;
        auto& left_neighbor = (*(--it_lower)).second;
        
        // If our point really is "between" the objective values of two points of that set
        
        if (new_point.obj_space[0] >= left_neighbor.obj_space[0] && 
            new_point.obj_space[0] <= right_neighbor.obj_space[0] && 
            new_point.obj_space[1] >= right_neighbor.obj_space[1] && 
            new_point.obj_space[1] <= left_neighbor.obj_space[1]) {
          
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