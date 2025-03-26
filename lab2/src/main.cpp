#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#define ITER 100

long double f(long double x) {
  return x * x * x / (6 * (std::sinh(x) - x));
}

int fac(int x) {
  return x < 2 ? 1 : x * fac(x-1);
}

long double f_tylor(long double x) {
  long double sum = 0.0L;
  long double term = x;

  for (int i = 1; i <= ITER; i++) {
    // z rozwiniecia w szereg tylora
    term *= x * x / fac(2 * i + 1);

    if (term < std::numeric_limits<long double>::epsilon())
      break;

    sum += term;
  }

  return (x * x * x) / (6.0L * sum);
}

int main() {
  long double logx, x, fx, e, et;
  std::string line;
  std::ifstream file("dane.txt");

  if (!file) {
    return 1;
  }

  while (std::getline(file, line)) {
    std::stringstream ss(line);

    ss >> logx >> x >> fx;

    e = std::fabsl((f(x) - fx) / fx);
    et = std::fabsl((f_tylor(x) - fx) / fx);

    //std::cout << logx << "\t" << std::log10(e) << "\t" << std::log10(et) << std::endl;
    std::cout << x << "\t" << fx << "\t" << f(x) << "\t" << f_tylor(x) << std::endl;
  }

  file.close();

  return 0;
}
