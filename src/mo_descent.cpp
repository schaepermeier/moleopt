#include "mo_descent.h"
#include "utils.h"
#include <cmath>

using namespace std;

double_vector mo_steepest_descent_direction(const std::vector<double_vector>& gradients) {
  int n_objectives = gradients.size();
  
  if (n_objectives != 2) {
    print("Cannot compute descent direction for n_objectives != 2");
    return double_vector();
  }
  
  // Steepest descent = Highest HV Improvement for given step size
  // Indicated by sqrt(HV)
  
  double_vector mog = -0.5 * (normalize(gradients[0]) + normalize(gradients[1]))
    * sqrt(norm(gradients[0])) * sqrt(norm(gradients[1]));
  
  return mog;
}

double_vector mo_steepest_descent_direction(const std::vector<double_vector>& gradients,
                                            const double_vector ref_point,
                                            const double_vector current_y) {
  double current_improvement = sqrt(max(0.0, ref_point[0] - current_y[0])) *
                               sqrt(max(0.0, ref_point[1] - current_y[1]));
  
  double_vector mog;
  
  if (current_improvement > 0) {
    // Hypervolume gradient direction:
    
    mog = (ref_point[1] - current_y[1]) * gradients[0] +
          (ref_point[0] - current_y[0]) * gradients[1];
    
    // Normalize and get descent direction
    mog = -mog / (2 * current_improvement);
  } else {
    mog = mo_steepest_descent_direction(gradients);
  }
  
  return mog;
}


