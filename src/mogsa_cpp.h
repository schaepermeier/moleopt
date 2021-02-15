#ifndef MOGSA_CPP_H
#define MOGSA_CPP_H

#include "vector_utils.h"
#include "types.h"
#include <map>

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
    double maximum_explore_set);

#endif
