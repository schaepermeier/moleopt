#ifndef SET_UTILS_H
#define SET_UTILS_H

#include "vector_utils.h"
#include "types.h"

#include <set>

int check_duplicated_set(const vector<efficient_set>&,
                         const evaluated_point&,
                         double epsilon);

vector<evaluated_point> sort_by_nondominated(vector<evaluated_point> point_set);

void insert_into_set(efficient_set& set,
                     const evaluated_point& new_point);

void refine_sets(vector<efficient_set>& sets,
                 double rel_hv_target,
                 const optim_fn& fn,
                 const corrector_fn& descent_fn,
                 set<double_vector>& nondominated_points);

bool is_nondominated(const set<double_vector>& nondominated_points,
                     const double_vector& obj_vector);

bool insert_nondominated(set<double_vector>& nondominated_points,
                         const double_vector& obj_vector);

bool point_dominated_by_set(const double_vector& a, const set<double_vector>& s);

#endif
