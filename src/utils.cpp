#include "utils.h"

// === Logging Helper ===

void print(std::string message) {
  std::cout << message << std::endl;
}

void print(double a) {
  std::cout << a << std::endl;
}

void print_vector(double_vector v) {
  for (const auto& el : v) std::cout << el << " ";
  
  std::cout << std::endl;
}

// === RNG ===

std::default_random_engine generator;
std::normal_distribution<double> distribution;

double random_double() {
  return distribution(generator);
}
