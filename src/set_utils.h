#ifndef SET_UTILS_H
#define SET_UTILS_H

#include "vector_utils.h"
#include "types.h"

#include <set>

int check_duplicated_set(const vector<efficient_set>&,
                         const evaluated_point&,
                         double);

void insert_into_set(efficient_set& set,
                     const evaluated_point& new_point);

void refine_sets(vector<efficient_set>& sets,
                 double rel_hv_target,
                 optim_fn fn,
                 corrector_fn descent_fn,
                 set<double_vector>& nondominated_points);

bool is_nondominated(set<double_vector>& nondominated_points, double_vector& obj_vector);

void insert_nondominated(set<double_vector>& nondominated_points, double_vector& obj_vector);

#endif
