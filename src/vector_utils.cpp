#include "vector_utils.h"
#include <cmath>
#include <algorithm>

double_vector operator+(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double_vector result;
  result.reserve(a.size());
  
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<double>());
  return result;
}

double_vector operator-(const double_vector& a, const double_vector& b) {
  assert(a.size() == b.size());
  
  double_vector result;
  result.reserve(a.size());
  
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<double>());
  return result;
}

double_vector operator*(const double_vector& a, double scalar) {
  double_vector result;
  result.reserve(a.size());
  
  std::transform(a.begin(), a.end(), std::back_inserter(result), [&](double v){return scalar * v;});
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
  
  std::transform(a.begin(), a.end(), std::back_inserter(result), [&](double v){return v / divisor;});
  return result;
}

double_vector operator/(double a, const double_vector& divisor) {
  assert(divisor != 0);
  
  double_vector result;
  result.reserve(divisor.size());
  
  std::transform(divisor.begin(), divisor.end(), std::back_inserter(result), [&](double v){return a / v;});
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

double_vector ensure_boundary(const double_vector& vector, const double_vector& lower, const double_vector& upper) {
  assert(vector.size() == lower.size());
  assert(vector.size() == upper.size());
  
  double_vector result(vector.size());
  
  for (int i = 0; i < vector.size(); i++) {
    result[i] = std::min(std::max(vector[i], lower[i]), upper[i]);
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
