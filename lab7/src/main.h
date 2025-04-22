#ifndef MAIN_H
#define MAIN_H

#include "matrix.h"

static const int MAX_ITERATIONS = 100;
static const long double TOLX = 1.0e-6l;
static const long double TOLF = 1.0e-6l;

int main();

void Jacobi(const Matrix& a, const Vector& b, const Vector& xStart);

void GaussSeidel(const Matrix& a, const Vector& b, const Vector& xStart);

void SOR(const Matrix& a, const Vector& b, const Vector& xStart,
    long double omega);

#endif
