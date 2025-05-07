#ifndef THOMAS_HPP
#define THOMAS_HPP

#include <vector>

typedef std::vector<double> Vector;

int main();

void CalculateEta(const Vector& u, const Vector& d, const Vector& l,
    Vector& eta);

void CalculateR(const Vector& eta, const Vector& l, const Vector& b,
    Vector& r);

void Solve(const Vector& u, const Vector& eta, const Vector& r, Vector& x);

void PrintVector(const Vector& x);

#endif
