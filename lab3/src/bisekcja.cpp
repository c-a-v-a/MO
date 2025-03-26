#include <cmath>
#include <iostream>

#include "bisekcja.hpp"

int main() {
  std::cout << "=== f1 ===" << std::endl;
  bisekcja(f1, -1, 1);
  std::cout << "=== f2 ===" << std::endl;
  bisekcja(f2, -1., 1.);
  return 0;
}

void bisekcja(double (*f)(double), double a, double b) {
  int i = 0;

  double x;

  double e;
  double residuum;

  if (b < a) {
    std::cout << "b musi byc wieksze od a" << std::endl;
    return;
  }

  if (f(a) * f(b) > 0) {
    std::cout << "a i b maja ten sam znak" << std::endl;
    return;
  }

  do {
    x = (a + b) / 2.;
    e = (b - a) / 2.;
    residuum = f(x);

    std::cout << "xn: " << x << "\t";
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

    if (f(x) * f(a) >= 0) a = x;
    else b = x;
    i++;
  } while (i < MAX_ITERATIONS);

  std::cout << "maksymalna liczba iteracji" << std::endl;
}

double f1(double x) {
  return std::tanh(x) + 2 * (x - 1);
}

double f2(double x) {
  return std::sinh(x) + x / 4 - 1;
}
