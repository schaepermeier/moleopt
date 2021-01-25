#ifndef MO_DESCENT_H
#define MO_DESCENT_H

#include "vector_utils.h"

double_vector mo_steepest_descent_direction(const std::vector<double_vector>& gradients);

double_vector mo_steepest_descent_direction(const std::vector<double_vector>& gradients,
                                            const double_vector ref_point,
                                            const double_vector current_y);

#endif
