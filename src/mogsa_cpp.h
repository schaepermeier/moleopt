#ifndef MOGSA_CPP_H
#define MOGSA_CPP_H

#include "vector_utils.h"
#include "types.h"
#include <map>

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
                const optim_fn& mo_function,
                const gradient_fn& gradient_function,
                const corrector_fn& descent_function,
                const explore_set_fn& explore_set_function,
                const vector<double_vector>& starting_points,
                const double_vector& lower_bounds,
                const double_vector& upper_bounds,
                double epsilon_explore_set,
                double epsilon_initial_step_size,
                double maximum_explore_set);

tuple<vector<efficient_set>, vector<tuple<int, int>>> run_mogsa(
                const optim_fn& mo_function,
                const vector<double_vector>& starting_points,
                const double_vector& lower,
                const double_vector& upper,
                double epsilon_gradient,
                double epsilon_explore_set,
                double epsilon_initial_step_size,
                double max_explore_set);

#endif
