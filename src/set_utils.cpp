#include "set_utils.h"

int check_duplicated_set(std::vector<efficient_set> local_sets,
                         evaluated_point current_point,
                         double eps_initial_step_size) {
  int containing_set = -1;
  
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
            return containing_set;
          }
        }
      }
    }
    
    for (auto& [f1_val, point] : set) {
      if (norm(current_point.dec_space - point.dec_space) < eps_initial_step_size) {
        return containing_set;
      }
    }
  }
  
  return -1;
}