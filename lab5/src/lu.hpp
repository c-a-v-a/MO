#ifndef LU_HPP
#define LU_HPP

#include <vector>

typedef std::vector<long double> Vector;
typedef std::vector<std::vector<long double>> Matrix;

bool IsSquare(const Matrix& x);

Vector GetColumn(const Matrix& x, int column);

int MaxValueIndex(const Vector& x, int start);

void LUFactorization(Matrix& l, Matrix& u, Vector& p);

void PrintVector(const Vector& x);

void PrintMatrix(const Matrix& x);

void PrintPVector(const Vector& x, const Vector& p);

void PrintPMatrix(const Matrix& x, const Vector& p);

void ZeroMatrix(Matrix& x, int size);

void ZeroVector(Vector& x, int size);

void PermutationVector(Vector& p, int size);

void VectorSubtract(Vector& x, const Vector& y, long double k = 1.0L);

void CalculateY(const Matrix& l, Vector& y, const Vector& b, const Vector& p);

void CalculateX(const Matrix& u, Vector& x, const Vector& y, const Vector& p);

void Solve(Matrix& a, Vector& x, const Vector& b);

#endif
