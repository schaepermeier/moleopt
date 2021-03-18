#include "utils.h"

// === Logging Helper ===

bool logging_enabled = true;

void print_info(const string& message) {
  cout << message << endl;
}

void print(const string& message) {
  if (logging_enabled) {
    cout << message << endl;
  }
}

void print(double a) {
  if (logging_enabled) {
    cout << a << endl;
  }
}

void print(const double_vector& v) {
  if (logging_enabled) {
    for (const auto& el : v) cout << el << " ";
    
    cout << endl;
  }
}

// === RNG ===

default_random_engine generator;
normal_distribution<double> distribution;

double random_double() {
  return distribution(generator);
}
