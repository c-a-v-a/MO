#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "thomas.h"

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
    if (x.at(i) == 0.0)
      x.at(i) = (r.at(i) - u.at(i) * x.at(i + 1)) / eta.at(i);
}

void PrintVector(const Vector& x) {
  for (auto y : x)
    std::cout << y << " ";

  std::cout << std::endl;
}
