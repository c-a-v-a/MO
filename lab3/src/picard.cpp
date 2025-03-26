#include <cmath>
#include <iostream>

#include "picard.hpp"

int main() {
  std::cout << "=== f1 ===" << std::endl;
  picard(f1, phi1);
  std::cout << "=== f2 ===" << std::endl;
  picard(f2, phi2);
  return 0;
}

void picard(double (*f)(double), double (*phi)(double)) {
  int i = 0;

  double x;
  double xn;

  double e;
  double residuum;

  x = (double)0;

  do {
    if (std::fabs(phi(x)) >= 1) {
      std::cout << "rozbieznosc" << std::endl;
      return;
    }

    xn = phi(x);
    e = xn - x;
    residuum = f(xn);

    std::cout << "xn: " << xn << "\t";
    std::cout << "e: " << e << "\t";
    std::cout << "res: " << residuum << std::endl;

    if (std::fabs(e) <= TOLX) {
      std::cout << "blad mniejszy od TOLX" << std::endl;
      return;
    }

    if (std::fabs(residuum) <= TOLF) {
      std::cout << "residuum mniejsze od TOLF" << std::endl;
      return;
    }

    i++;
    x = xn;
  } while (i < MAX_ITERATIONS);

  std::cout << "osiagnieto maksymalna liczbe iteracji" << std::endl;
}

double phi1(double x) {
  return 0.5 - std::pow(std::tanh(x), 2);
}

double f1(double x) {
  return std::tanh(x) + 2 * (x - 1);
}

double phi2(double x) {
  return std::asinh(-1 * x / 4 + 1);
}

double f2(double x) {
  return std::sinh(x) + x / 4 - 1;
}
