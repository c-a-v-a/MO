#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#define ITER 250

double f(double x) {
  return x * x * x / (6 * (std::sinh(x) - x));
}

int fac(int x) {
  return x <= 2 ? 1 : x * fac(x-1);
}

double f_tylor(double x, int iter) {
  double result = x * x * x;
  double sum;

  for (int i = 0; i < iter; i++) {
    // z rozwiniecia w szereg tylora
    sum += std::pow(x, 2 * i + 1) / fac(2 * i + 1);
  }

  sum -= x;
  sum *= 6;

  result /= sum;

  return result;
}

int main() {
  double logx, x, fx, e;
  std::string line;
  std::ifstream file("dane.txt");

  if (!file) {
    return 1;
  }

  while (std::getline(file, line)) {
    std::stringstream ss(line);

    ss >> logx >> x >> fx;

    e = std::fabs((f(x) - fx) / fx);

    std::cout << logx << "\t" << std::log10(e) << std::endl;
  }

  file.close();

  return 0;
}
