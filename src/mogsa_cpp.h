#ifndef MOGSA_CPP_H
#define MOGSA_CPP_H

#include "vector_utils.h"
#include <map>

std::tuple<evaluated_point, std::vector<evaluated_point>> explore_efficient_set(evaluated_point current_point, int objective);

std::tuple<std::vector<std::map<double, evaluated_point>>,
           std::vector<std::tuple<int, int>>> run_mogsa(optim_fn f, std::vector<double_vector> starting_points, double_vector lower, double_vector upper,
               double epsilon_gradient, double epsilon_explore_set, double epsilon_initial_step_size, double maximum_explore_set);

#endif
