#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "thomas.hpp"

int main() {
  Vector u = { -1.0L, -3.0L, 5.0L, -7.0L };
  Vector d = { 100.0L, 200.0L, 300.0L, 200.0L, 100.0L };
  Vector l = { 2.0L, 4.0L, -6.0L, -8.0L };
  Vector b = { 199.0L, 195.0L, 929.0L, 954.0L, 360.0L };
  Vector x;
  Vector eta;
  Vector r;

  std::cout << "=== VECTOR U ===" << std::endl;
  PrintVector(u);
  std::cout << std::endl;

  std::cout << "=== VECTOR D ===" << std::endl;
  PrintVector(d);
  std::cout << std::endl;

  std::cout << "=== VECTOR L ===" << std::endl;
  PrintVector(l);
  std::cout << std::endl;

  CalculateEta(u, d, l, eta);

  std::cout << "=== VECTOR ETA ===" << std::endl;
  PrintVector(eta);
  std::cout << std::endl;

  CalculateR(eta, l, b, r);

  std::cout << "=== VECTOR R ===" << std::endl;
  PrintVector(r);
  std::cout << std::endl;

  Solve(u, eta, r, x);

  std::cout << "=== VECTOR X ===" << std::endl;
  PrintVector(x);
  std::cout << std::endl;
}

void CalculateEta(const Vector& u, const Vector& d, const Vector& l,
    Vector& eta) {
  if (u.size() != l.size() || d.size() - 1 != u.size())
    throw std::runtime_error("Invalid vector sizes");

  eta.resize(d.size());

  eta.at(0) = d.at(0);

  for (int i = 1; i < d.size(); i++)
    eta.at(i) = d.at(i) - l.at(i - 1) * u.at(i - 1) / eta.at(i - 1);
}

void CalculateR(const Vector& eta, const Vector& l, const Vector& b,
    Vector& r) {
  if (eta.size() != b.size() || eta.size() - 1 != l.size())
    throw std::runtime_error("Invalid vector sizes"); 

  r.resize(b.size());

  r.at(0) = b.at(0);

  for (int i = 1; i < b.size(); i++)
    r.at(i) = b.at(i) - l.at(i - 1) * r.at(i - 1) / eta.at(i - 1);
}

void Solve(const Vector& u, const Vector& eta, const Vector& r,
    Vector& x) {
  if (eta.size() != r.size() || eta.size() - 1 != u.size())
    throw std::runtime_error("Invalid vector sizes"); 

  x.resize(eta.size());

  x.back() = r.back() / eta.back();

  for (int i = eta.size() - 2; i >= 0; i--)
    x.at(i) = (r.at(i) - u.at(i) * x.at(i + 1)) / eta.at(i);
}

void PrintVector(const Vector& x) {
  for (auto y : x)
    std::cout << y << " ";

  std::cout << std::endl;
}
