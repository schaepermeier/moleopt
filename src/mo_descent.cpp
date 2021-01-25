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
                                            const double_vector& ref_point,
                                            const double_vector& current_y) {
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

gradient_fn create_gradient_fn(optim_fn fn,
                               const double_vector& lower,
                               const double_vector& upper,
                               const std::string& method,
                               double eps_gradient) {
  
  gradient_fn g = [fn, lower, upper, method, eps_gradient](const evaluated_point& point) {
    std::vector<double_vector> gradients(2);
    int d = point.dec_space.size();
    
    // Initialize gradients
    
    gradients[0] = double_vector(d, 0);
    gradients[1] = double_vector(d, 0);
    
    for (int iter_d = 0; iter_d < d; iter_d++) {
      double_vector lower_d(point.dec_space);
      double_vector upper_d(point.dec_space);
      double_vector fn_lower;
      double_vector fn_upper;
      
      if (method == "twosided") {
        lower_d[iter_d] = lower_d[iter_d] - eps_gradient;
        upper_d[iter_d] = upper_d[iter_d] + eps_gradient;
      } else {
        if (point.dec_space[iter_d] == upper[iter_d]) {
          // point == upper
          upper_d[iter_d] = point.dec_space[iter_d];
          lower_d[iter_d] = point.dec_space[iter_d] - eps_gradient;
        } else if (point.dec_space[iter_d] == lower[iter_d]){
          // point == lower
          upper_d[iter_d] = point.dec_space[iter_d] + eps_gradient;
          lower_d[iter_d] = point.dec_space[iter_d];
        } else {
          // Randomized forward or backward differences
          if (random_double() >= 0) {
            upper_d[iter_d] = point.dec_space[iter_d];
            lower_d[iter_d] = point.dec_space[iter_d] - eps_gradient;
          } else {
            upper_d[iter_d] = point.dec_space[iter_d] + eps_gradient;
            lower_d[iter_d] = point.dec_space[iter_d];
          }
        }
      }
      
      lower_d = ensure_boundary(lower_d, lower, upper);
      upper_d = ensure_boundary(upper_d, lower, upper);
      
      if (lower_d == point.dec_space) {
        fn_lower = point.obj_space;
      } else {
        fn_lower = fn(lower_d);
      }
      
      if (upper_d == point.dec_space) {
        fn_upper = point.obj_space;
      } else {
        fn_upper = fn(upper_d);
      }
      
      double length = norm(lower_d - upper_d);
      
      double_vector gradient_components = (fn_upper - fn_lower) / length;
      
      gradients[0][iter_d] = gradient_components[0];
      gradients[1][iter_d] = gradient_components[1];
    }
    
    return gradients;
  };
  
  return g;
}
