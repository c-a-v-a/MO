#include <cmath>
#include <iostream>
#include <stdexcept>

#include "matrix.h"

void CreateVector(Vector& x, int size) {
  x.size = size;
  x.vector = new long double[size];
}

void CreateMatrix(Matrix& x, int size) {
  x.size = size;
  x.matrix = new long double*[size];
  
  for (int i = 0; i < size; i++)
    x.matrix[i] = new long double[size];
}

void FreeVector(Vector& x) {
  delete[] x.vector;
}

void FreeMatrix(Matrix& x) {
  for (int i = 0; i < x.size; i++)
    delete[] x.matrix[i];

  delete[] x.matrix;
}

void IdentityMatrix(Matrix& x) {
  for (int i = 0; i < x.size; i++)
    for (int j = 0; j < x.size; j++)
      x.matrix[i][j] = i == j ? 1.0L : 0.0L;
}

void LUSplit(Matrix& ld, Matrix& u) {
  for (int i = 0; i < u.size; i++)
    for (int j = 0; j < u.size; j++)
      u.matrix[i][j] = 0.0L;

  for (int i = 0; i < ld.size; i++) {
    for (int j = i + 1; j < ld.size; j++) {
      u.matrix[i][j] = ld.matrix[i][j];
      ld.matrix[i][j] = 0.0L;
    }
  }
}

void PrintVector(const Vector& x) {
  for (int i = 0; i < x.size; i++)
    std::cout << x.vector[i] << " ";

  std::cout << std::endl;
}

void PrintMatrix(const Matrix& x) {
  for (int i = 0; i < x.size; i++) {
    for (int j = 0; j < x.size; j++)
      std::cout << x.matrix[i][j] << " ";

    std::cout << std::endl;
  }
}

void DiagonalInverse(const Matrix& x, Matrix& y) {
  for (int i = 0; i < x.size; i++) {
    if (std::fpclassify(x.matrix[i][i]) == FP_ZERO)
      throw std::runtime_error("Matrix is not invertible");

    y.matrix[i][i] = 1.0L / x.matrix[i][i];
  }
}

void LowerTriangularInverse(const Matrix& x, Matrix& y) {
  IdentityMatrix(y); 

  for (int i = 0; i < x.size; i++) {
    y.matrix[0][i] = y.matrix[0][i] / x.matrix[0][0];

    for (int j = 1; j < x.size; j++) {
      long double sum = y.matrix[j][i];

      for (int k = 0; k < j; k++)
        sum -= x.matrix[j][k] * y.matrix[k][i];

      y.matrix[j][i] = sum / x.matrix[j][j];
    }
  }
}

void MVMultiplication(const Matrix& x, const Vector& y, Vector& z) {
  for (int i = 0; i < x.size; i++) {
    z.vector[i] = 0.0L;

    for (int j = 0; j < x.size; j++)
      z.vector[i] += x.matrix[i][j] * y.vector[j];
  }
}

void MMMultiplication(const Matrix& x, const Matrix& y, Matrix& z) {
  for (int i = 0; i < x.size; i++) {
    for (int j = 0; j < x.size; j++) {
      z.matrix[i][j] = 0.0L;

      for (int k = 0; k < x.size; k++)
        z.matrix[i][j] += x.matrix[i][k] * y.matrix[k][j];
    }
  }
}

void MSMultiplication(Matrix& x, long double s) {
  for (int i = 0; i < x.size; i++)
    for (int j = 0; j < x.size; j++)
      x.matrix[i][j] = s * x.matrix[i][j];
}

void VectorAdd(Vector& x, const Vector& y) {
  for (int i = 0; i < x.size; i++)
    x.vector[i] += y.vector[i];
}

void VectorSubtract(Vector& x, const Vector& y) {
  for (int i = 0; i < x.size; i++)
    x.vector[i] = y.vector[i] - x.vector[i];
}

long double MatrixMaxNorm(const Matrix& x) {
  long double norm = 0.0L;

  for (int i = 0; i < x.size; i++) {
    long double sum = 0.0L;

    for (int j = 0; j < x.size; j++)
      sum += std::fabs(x.matrix[i][j]);

    if (sum > norm)
      norm = sum;
  }

  return norm;
}

long double VectorMaxNorm(const Vector& x) {
  long double norm = 0.0L;

  for (int i = 0; i < x.size; i++)
    norm = x.vector[i] > norm ? x.vector[i] : norm;

  return norm;
}
