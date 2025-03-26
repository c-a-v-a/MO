#include <cmath>
#include <iostream>

#include "sieczne.hpp"

int main() {
  std::cout << "=== f1 ===" << std::endl;
  sieczne(f1, -1, 0);
  std::cout << "=== f2 ===" << std::endl;
  sieczne(f2, -1, 0);
  return 0;
}

void sieczne(double (*f)(double), double x, double xn) {
  int i = 0;

  double xn2;

  double e;
  double residuum;

  do {
    xn2 = xn - f(xn) / ((f(xn) - f(x)) / (xn - x));
    e = xn2 - xn;
    residuum = f(xn2);

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
    xn = xn2;
  } while (i < MAX_ITERATIONS);

  std::cout << "maksymalna liczba iteracji" << std::endl;
}

double f1(double x) {
  return std::tanh(x) + 2 * (x - 1);
}

double f2(double x) {
  return std::sinh(x) + x / 4 - 1;
}

