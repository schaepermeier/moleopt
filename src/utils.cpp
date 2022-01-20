#include "utils.h"

// === Logging Helper ===

MOLE_LOG_LEVEL_ENUM MOLE_LOG_LEVEL = MOLE_LOG_LEVEL_INFO;

void print_info(const string& message) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_INFO) {
    cout << message << endl;
  }
}

void print_info(double a) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_INFO) {
    cout << a << endl;
  }
}

void print_info(const double_vector& v) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_INFO) {
    for (const auto& el : v) cout << el << " ";
    
    cout << endl;
  }
}

void print(const string& message) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_DEBUG) {
    cout << message << endl;
  }
}

void print(double a) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_DEBUG) {
    cout << a << endl;
  }
}

void print(const double_vector& v) {
  if (MOLE_LOG_LEVEL >= MOLE_LOG_LEVEL_DEBUG) {
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
