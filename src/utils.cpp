#include "utils.h"

// === Logging Helper ===

void print(string message) {
  cout << message << endl;
}

void print(double a) {
  cout << a << endl;
}

void print(double_vector v) {
  for (const auto& el : v) cout << el << " ";
  
  cout << endl;
}

// === RNG ===

default_random_engine generator;
normal_distribution<double> distribution;

double random_double() {
  return distribution(generator);
}
