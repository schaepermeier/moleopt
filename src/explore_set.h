#ifndef EXPLORE_SET_H
#define EXPLORE_SET_H

#include "types.h"

explore_set_fn get_explore_set_fn(
        const optim_fn& fn,
        const gradient_fn& grad_fn,
        const corrector_fn& descent_fn,
        const double_vector& lower,
        const double_vector& upper,
        double explore_step_min,
        double explore_step_max,
        double explore_angle_max,
        double explore_scale_factor);

#endif
