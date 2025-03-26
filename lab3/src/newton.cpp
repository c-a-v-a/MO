#include <cmath>
#include <iostream>

#include "newton.hpp"

int main() {
  std::cout << "=== f1 ===" << std::endl;
  newton(f1, df1);
  std::cout << "=== f2 ===" << std::endl;
  newton(f2, df2);
  return 0;
}

void newton(double (*f)(double), double (*df)(double)) {
  int i = 0;

  double x;
  double xn;

  double e;
  double residuum;

  x = 0.7;

  do {
    xn = x - (f(x) / df(x));
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

  std::cout << "maksymalna liczba iteracji" << std::endl;
}

double f1(double x) {
  return std::tanh(x) + 2 * (x - 1);
}

double df1(double x) {
  return std::tanh(x) * std::tanh(x) + 1;
}

double f2(double x) {
  return std::sinh(x) + x / 4 - 1;
}

double df2(double x) {
  return std::cosh(x) + 0.25;
}
