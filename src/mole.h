#ifndef MOLE_H
#define MOLE_H

#include "vector_utils.h"
#include "types.h"
#include <map>

tuple<vector<efficient_set>,
      vector<tuple<int, int>>,
      vector<long>> run_mole(
                const optim_fn& mo_function,
                const gradient_fn& gradient_function,
                const corrector_fn& descent_function,
                const explore_set_fn& explore_set_function,
                const vector<double_vector>& starting_points,
                const double_vector& lower_bounds,
                const double_vector& upper_bounds,
                int max_local_sets,
                double explore_step_min,
                int refine_after_nstarts,
                double refine_hv_target,
                long max_budget,
                long* used_budget);

tuple<vector<efficient_set>,
      vector<tuple<int, int>>,
      vector<long>> run_mole(
                const optim_fn& mo_function,
                const vector<double_vector>& starting_points,
                const double_vector& lower,
                const double_vector& upper,
                int max_local_sets,
                double epsilon_gradient,
                double descent_direction_min,
                double descent_step_min,
                double descent_step_max,
                double descent_scale_factor,
                double descent_armijo_factor,
                int descent_history_size,
                int descent_max_iter,
                double explore_step_min,
                double explore_step_max,
                double explore_angle_max,
                double explore_scale_factor,
                int refine_after_nstarts,
                double refine_hv_target);

#endif
