#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
  long double* vector;
  int size;
} Vector;

typedef struct {
  long double** matrix;
  int size;
} Matrix;

void CreateVector(Vector& x, int size);

void CreateMatrix(Matrix& x, int size);

void FreeVector(Vector& x);

void FreeMatrix(Matrix& x);

void IdentityMatrix(Matrix& x);

void LUSplit(Matrix& ld, Matrix& u);

void PrintVector(const Vector& x);

void PrintMatrix(const Matrix& x);

void DiagonalInverse(const Matrix& x, Matrix& y);

void LowerTriangularInverse(const Matrix& x, Matrix& y);

void MVMultiplication(const Matrix& x, const Vector& y, Vector& z);

void MMMultiplication(const Matrix& x, const Matrix& y, Matrix& z);

void MSMultiplication(Matrix& x, long double s);

void VectorAdd(Vector& x, const Vector& y);

void VectorSubtract(Vector& x, const Vector& y);

long double MatrixMaxNorm(const Matrix& x);

long double VectorMaxNorm(const Vector& x);

#endif
