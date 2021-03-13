#ifndef SET_UTILS_H
#define SET_UTILS_H

#include "vector_utils.h"
#include "types.h"

int check_duplicated_set(const vector<efficient_set>&,
                         const evaluated_point&,
                         double);

void insert_into_set(efficient_set& set,
                     const evaluated_point& new_point);

void refine_sets(vector<efficient_set>& sets,
                 double rel_hv_target,
                 optim_fn fn,
                 corrector_fn descent_fn);

#endif
