#include <cmath>
#include <fstream>
#include <iostream>

#include "thomas.h"

#include "threepoint.h"

int main() {
  std::ofstream errfile("threepoint_error.csv");
  double N = 1.5e8;

  for (double i = 10.0; i < N; i *= 1.5) {
    double h = (b - a) / i;
    double error = -1.0;
    Vector y;

    threepoint(i, y);

    for (int j = 0; j < y.size(); j++) {
      y.at(j) = std::fabs(y.at(j) - u(a + j * h));

      if (y.at(j) > error)
        error = y.at(j);
    }

    std::cout << std::log10(error) << " " << std::log10(h) << std::endl;
    errfile << std::log10(error) << " " << std::log10(h) << std::endl;
  }

  errfile.close();

  std::ofstream outfile("threepoint.csv");
  double h = (b - a) / 1.0e6;
  Vector y;

  threepoint(1.0e6, y);

  for (int j = 0; j < y.size(); j++)
    outfile << (a + j * h) << " " << y.at(j) << std::endl;

  outfile.close();
}

void threepoint(double N, Vector& y) {
  double h = (b - a) / N;

  // vectors for thomas algorithm
  Vector u;
  Vector d;
  Vector l;
  Vector bv;
  Vector eta;
  Vector rv;

  // the middle part of matrix is constant, since p,q,r and h are constant.
  const double l_term = p / (h * h) - q / (2 * h); 
  const double d_term = r - 2 * p / (h * h);
  const double u_term = p / (h * h) + q / (2 * h);

  d.push_back(beta - alpha / h);
  u.push_back(alpha / h);
  bv.push_back(-1.0 * gama);

  for (int i = 1; i < N; i++) {
    l.push_back(l_term);
    d.push_back(d_term);
    u.push_back(u_term);
    bv.push_back(-1 * s(a + i * h));
  }

  l.push_back(-1.0 * phi / h);
  d.push_back(phi / h + psi);
  bv.push_back(-1 * theta);

  // thomas algorithm
  CalculateEta(u, d, l, eta);
  CalculateR(eta, l, bv, rv);
  Solve(u, eta, rv, y);
}
