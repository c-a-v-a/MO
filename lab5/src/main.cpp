#include "lu.hpp"

int main() {
  Matrix a = {
    { 5.0L, 4.0L, 3.0L, 2.0L, 1.0L },
    { 10.0L, 8.0L, 7.0L, 6.0L, 5.0L },
    { -1.0L, 2.0L, -3.0L, 4.0L, -5.0L },
    { 6.0L, 5.0L, -4.0L, 3.0L, -2.0L },
    { 1.0L, 2.0L, 3.0L, 4.0L, 5.0L }
  };
  Vector x;
  Vector b = { 37.0L, 99.0L, -9.0L, 12.0L, 53.0L };

  Solve(a, x, b);
}
