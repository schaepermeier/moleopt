#ifndef EXPLORE_SET_H
#define EXPLORE_SET_H

#include "types.h"

explore_set_fn get_explore_set_fn(
    const optim_fn& fn,
    const gradient_fn& grad_fn,
    const corrector_fn& descent_fn,
    const double_vector& lower,
    const double_vector& upper,
    double eps_explore_set,
    double max_explore_set);

#endif
