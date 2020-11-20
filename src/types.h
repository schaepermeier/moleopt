#ifndef TYPES_H
#define TYPES_H

#include <vector>

typedef std::vector<double> double_vector;
typedef std::function<std::vector<double>(std::vector<double>)> optim_fn;

struct evaluated_point {
  double_vector dec_space;
  double_vector obj_space;
  std::vector<double_vector> gradients;
  
  // bool operator<(const evaluated_point& rhs) const {
  //   int i = 0;
  //   int max_index = std::max(this->obj_space.size(), rhs.obj_space.size()) - 1;
  //   
  //   if (this->obj_space[i] == rhs.obj_space[i]) {
  //     if (i < max_index) i++;
  //   }
  //   
  //   return this->obj_space[i] < rhs.obj_space[i];
  // }
};

#endif