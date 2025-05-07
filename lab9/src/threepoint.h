#ifndef THREEPOINT_H
#define THREEPOINT_H

#include "thomas.h"

// interval
static const double a = 0.0;
static const double b = 1.0;

// given functions
static const double p = 1.0;
static const double q = 2.0;
static const double r = 4.0;
double s(double x) {
  return x * x * x / 2;
}

// boundry conditions
// for 'a' boundry
static const double alpha = 0.0;
static const double beta = 1.0;
static const double gama = -2.0;

// for 'b' boundry
static const double phi = 0.0;
static const double psi = 1.0;
static const double theta = 2.0;

void threepoint(double N, Vector& y);

double u(double x) {
    const double sqrt5 = std::sqrt(5.0);
    const double coth_sqrt5 = std::cosh(sqrt5) / std::sinh(sqrt5);

    double term1 = 9.0;
    double term2 = -95.0 * std::exp((-1.0 - sqrt5) * (-1.0 + x));
    double term3 = 55.0 * std::exp((-1.0 + sqrt5) * x);
    double term4 = 95.0 * std::exp(1.0 + sqrt5 + (-1.0 + sqrt5) * x);
    double term5 = -55.0 * std::exp(2.0 * sqrt5 - (1.0 + sqrt5) * x);
    double term6 = 2.0 * x * (6.0 + x * (3.0 + 2.0 * x));
    double term7 = -std::exp(2.0 * sqrt5) * (9.0 + 2.0 * x * (6.0 + x * (3.0 + 2.0 * x)));

    double numerator = (term1 + term2 + term3 + term4 + term5 + term6 + term7) * (-1.0 + coth_sqrt5);
    double result = -numerator / 64.0;

    return result;
}

#endif
