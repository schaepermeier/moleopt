#ifndef MO_DESCENT_H
#define MO_DESCENT_H

#include "vector_utils.h"
#include "utils.h"

gradient_fn create_gradient_fn(const optim_fn& fn,
                               const double_vector& lower,
                               const double_vector& upper,
                               const string& method,
                               double eps_gradient);

double_vector mo_steepest_descent_direction(const vector<double_vector>& gradients);

double_vector mo_steepest_descent_direction(const vector<double_vector>& gradients,
                                            const double_vector& ref_point,
                                            const double_vector& current_y);

corrector_fn create_two_point_stepsize_descent(const optim_fn& fn,
                                               const gradient_fn& grad_fn,
                                               double descent_direction_min,
                                               double descent_step_min,
                                               double descent_step_max,
                                               double descent_scale_factor,
                                               double descent_armijo_factor,
                                               int descent_history_size,
                                               int descent_max_iter,
                                               const double_vector& lower,
                                               const double_vector& upper);

#endif
