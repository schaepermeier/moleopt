#include "vector_utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

double_vector operator+(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), b.begin(), back_inserter(result), plus<double>());
  return result;
}

double_vector operator+(const double_vector& a, double scalar) {
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), back_inserter(result), [&](double v){return v + scalar;});
  return result;
}

double_vector operator-(const double_vector& a, double scalar) {
  return(a + (-scalar));
}

double_vector operator-(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), b.begin(), back_inserter(result), minus<double>());
  return result;
}

double_vector operator*(const double_vector& a, double scalar) {
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), back_inserter(result), [&](double v){return scalar * v;});
  return result;
}

double_vector operator*(double scalar, const double_vector& a) {
  return operator*(a, scalar);
}

double_vector operator-(const double_vector& a) {
  return operator*(a, -1.0);
}

double_vector operator/(const double_vector& a, double divisor) {
  assert(divisor != 0);
  
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), back_inserter(result), [&](double v){return v / divisor;});
  return result;
}

bool dominates(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  bool a_not_worse = true;
  bool a_strictly_better = false;
  
  for (int i = 0; i < a.size(); i++) {
    a_not_worse = a_not_worse && (a[i] <= b[i]);
    a_strictly_better = a_strictly_better || (a[i] < b[i]);
  }
  
  return a_not_worse && a_strictly_better;
}

bool strictly_dominates(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  bool a_strictly_better = true;
  
  for (int i = 0; i < a.size(); i++) {
    a_strictly_better = a_strictly_better && (a[i] < b[i]);
  }
  
  return a_strictly_better;
}

double norm(const double_vector& vector) {
  return sqrt(dot(vector, vector));
}
double square_norm(const double_vector& vector) {
  return dot(vector, vector);
}

double dot(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  return inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

double angle(const double_vector& a, const double_vector& b) {
  double enumerator = dot(a, b);
  double denominator = norm(a) * norm(b);
  
  if (denominator == 0) {
    return 180;
  }
  
  double rad = acos(enumerator / denominator);
  double angle = 180 * rad / M_PI;
  
  angle = max(min(angle, 180.0), 0.0);
  
  return angle;
}

void ensure_boundary(double_vector& vector, const double_vector& lower, const double_vector& upper) {
  assert(vector.size() == lower.size());
  assert(vector.size() == upper.size());

  for (int i = 0; i < vector.size(); i++) {
    vector[i] = min(max(vector[i], lower[i]), upper[i]);
  }
}

double_vector normalize(const double_vector& vector) {
  double length = norm(vector);
  
  if (length == 0.0) {
    return 0.0 * vector;
  } else {
    return vector / length;
  }
}

double compute_improvement(const double_vector& obj_space, const double_vector& ref_point) {
  
  if (ref_point.size() != 2) {
    throw "Improvement only supported for dim == 2";
  }
  
  if (dominates(obj_space, ref_point)) {
    return sqrt(max(0.0, ref_point[0] - obj_space[0])) *
           sqrt(max(0.0, ref_point[1] - obj_space[1]));
  } else if (obj_space == ref_point) {
    return 0;
  } else {
    return -inf;
  }
}

void project_feasible_direction(double_vector& search_direction, const double_vector& current_position,
                                const double_vector& lower, const double_vector& upper) {
  assert(search_direction.size() == lower.size());
  assert(search_direction.size() == upper.size());
  assert(search_direction.size() == current_position.size());

  for (int i = 0; i < search_direction.size(); i++) {
    if (current_position[i] == upper[i]) {
      // only negative direction allowed
      search_direction[i] = min(search_direction[i], 0.0);
    } else if (current_position[i] == lower[i]) {
      // only positive direction allowed
      search_direction[i] = max(search_direction[i], 0.0);
    } else {
      search_direction[i] = search_direction[i];
    }
  }
}
