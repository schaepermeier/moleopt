#ifndef EXPLORE_SET_H
#define EXPLORE_SET_H

#include "types.h"

explore_set_fn get_explore_set_fn(
    optim_fn fn,
    gradient_fn grad_fn,
    corrector_fn descent_fn,
    double_vector lower,
    double_vector upper,
    double eps_explore_set,
    double max_explore_set);

#endif
