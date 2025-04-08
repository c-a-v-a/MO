#include <cmath>
#include <iostream> 

#include "lu.hpp"

int main() {
  for (int i = -5; i > -500; i--) {
    long double e = std::powl(10.0L, i);

    Matrix A = {
      { 1.0L + e, 1.0L, 1.0L, 1.0L },
      { 1.0L, 1.0L + e, 1.0L, 1.0L },
      { 1.0L, 1.0L, 1.0L + e, 1.0L },
      { 1.0L, 1.0L, 1.0L, 1.0L + e }
    };
    Vector x;
    Vector b = { 6.0L + e, 6.0L + 2.0L * e, 6.0L + 2.0L * e, 6.0L + e };

    std::cout.precision(std::numeric_limits<long double>::max_digits10 - 1);
    std::cout << "=== wynik dla e=" << e << "===" << std::endl;

    Solve(A, x, b);

    if (!(x.at(0) == x.at(0))) break;
  }
}
