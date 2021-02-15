#include "vector_utils.h"
#include <cmath>
#include <algorithm>

using namespace std;

double_vector operator+(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double_vector result;
  result.reserve(a.size());
  
  transform(a.begin(), a.end(), b.begin(), back_inserter(result), plus<double>());
  return result;
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

double_vector operator/(double a, const double_vector& divisor) {
  assert(divisor != 0);
  
  double_vector result;
  result.reserve(divisor.size());
  
  transform(divisor.begin(), divisor.end(), back_inserter(result), [&](double v){return a / v;});
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

double square_norm(const double_vector& vector) {
  double square_norm = 0;
  
  for (double value : vector) {
    square_norm += value * value;
  }
  
  return square_norm;
}

double norm(const double_vector& vector) {
  return sqrt(square_norm(vector));
}

double dot(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double result = 0;
  
  for (int i = 0; i < a.size(); i++) {
    result += a[i] * b[i];
  }
  
  return result;
}

double angle(const double_vector& a, const double_vector& b) {
  double enumerator = dot(a, b);
  double denominator = norm(a) * norm(b);
  
  if (denominator == 0) {
    return 180;
  }
  
  double rad = acos(enumerator / denominator);
  double pi = atan(1) * 4;
  
  double angle = 180 * rad / pi;
  
  angle = max(min(angle, 180.0), 0.0);
  
  return angle;
}

double_vector ensure_boundary(const double_vector& vector, const double_vector& lower, const double_vector& upper) {
  assert(vector.size() == lower.size());
  assert(vector.size() == upper.size());
  
  double_vector result(vector.size());
  
  for (int i = 0; i < vector.size(); i++) {
    result[i] = min(max(vector[i], lower[i]), upper[i]);
  }
  
  return result;
}

double_vector normalize(const double_vector& vector) {
  double length = norm(vector);
  
  if (length == 0) {
    return vector;
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
