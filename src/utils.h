#ifndef UTILS_H
#define UTILS_H

#include "types.h"
#include <iostream>
#include <random>
#include <string>

// === Logging Helper ===

enum MOLE_LOG_LEVEL_ENUM {
  MOLE_LOG_LEVEL_NONE = 0,
  MOLE_LOG_LEVEL_INFO,
  MOLE_LOG_LEVEL_DEBUG
};

extern MOLE_LOG_LEVEL_ENUM MOLE_LOG_LEVEL;

void print_info(const string& message);

void print(const string& message);

void print(double a);

void print(const double_vector& v);

// === RNG ===

double random_double();

// === Infinity ===

double const inf = numeric_limits<double>::infinity();

#endif
