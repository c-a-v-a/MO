#include <cmath>
#include <iostream>

#include "main.h"
#include "matrix.h"

int main() {
  long double adata[5][5] = {
    { 50.0L, 5.0L, 4.0L, 3.0L, 2.0L },
    { 1.0L, 40.0L, 1.0L, 2.0L, 3.0L },
    { 4.0L, 5.0L, 30.0L, -5.0L, -4.0L },
    { -3.0L, -2.0L, -1.0L, 20.0L, 0.0L },
    { 1.0L, 2.0L, 3.0L, 4.0L, 30.0L}
  };
  long double bdata[5] = { 140.0L, 67.0L, 62.0L, 89.0L, 153.0L };
  long double xdata[5] = { 6.0L, 6.0L, 6.0L, 6.0L, 6.0L };
  long double omega = 0.5L;

  Matrix a;
  Vector b;
  Vector xStart;

  CreateMatrix(a, 5);
  CreateVector(b, 5);
  CreateVector(xStart, 5);

  for (int i = 0; i < a.size; i++) {
    for (int j = 0; j < a.size; j++)
      a.matrix[i][j] = adata[i][j];

    b.vector[i] = bdata[i];
    xStart.vector[i] = xdata[i];
  }

  std::cout << "=== JACOBI ===" << std::endl;
  Jacobi(a, b, xStart);
  std::cout << std::endl;

  std::cout << "=== GAUSS-SEIDEL ===" << std::endl;
  GaussSeidel(a, b, xStart);
  std::cout << std::endl;

  std::cout << "=== SOR ===" << std::endl;
  SOR(a, b, xStart, omega);
  std::cout << std::endl;

  FreeMatrix(a);
  FreeVector(b);
  FreeVector(xStart);
}

void Jacobi(const Matrix& a, const Vector& b, const Vector& xStart) {
  Vector x;
  Vector xn;
  Vector c;
  Vector r;
  Vector e;
  Matrix m;
  Matrix dInv;
  Matrix lu;

  int iter = 1;
  long double* temp;

  CreateVector(x, xStart.size);
  CreateVector(xn, xStart.size);
  CreateVector(c, xStart.size);
  CreateVector(r, xStart.size);
  CreateVector(e, xStart.size);
  CreateMatrix(m, a.size);
  CreateMatrix(dInv, a.size);
  CreateMatrix(lu, a.size);

  IdentityMatrix(dInv);
  
  DiagonalInverse(a, dInv);
  MVMultiplication(dInv, b, c);

  for (int i = 0; i < xStart.size; i++)
    x.vector[i] = xStart.vector[i];

  // L + U
  for (int i = 0; i < a.size; i++)
    for (int j = 0; j < a.size; j++)
      lu.matrix[i][j] = i != j ? a.matrix[i][j] : 0.0L;

  MMMultiplication(dInv, lu, m);
  MSMultiplication(m, -1.0L);

  if (MatrixMaxNorm(m) >= 1.0L) {
    std::cout << "Method not converging" << std::endl;

    FreeVector(x);
    FreeVector(xn);
    FreeVector(c);
    FreeVector(r);
    FreeVector(e);
    FreeMatrix(m);
    FreeMatrix(dInv);
    FreeMatrix(lu);

    return;
  }

  do {
    // calculate xn
    MVMultiplication(m, x, xn);
    VectorAdd(xn, c);

    // residuum
    MVMultiplication(a, xn, r);
    VectorSubtract(r, b);

    // err
    for (int i = 0; i < e.size; i++)
      e.vector[i] = std::fabs(xn.vector[i] - x.vector[i]);

    std::cout << "--- x" << iter << " ---" << std::endl;
    PrintVector(xn);

    std::cout << "--- residuum ---" << std::endl;
    PrintVector(r);
    std::cout << "MaxNorm: " << VectorMaxNorm(r) << std::endl;

    std::cout << "--- error estimator ---" << std::endl;
    PrintVector(e);
    std::cout << "MaxNorm: " << VectorMaxNorm(e) << std::endl;
    std::cout << std::endl;

    if (VectorMaxNorm(e) < TOLX && VectorMaxNorm(r) < TOLF) {
      std::cout << "Precision reached" << std::endl;
      break;
    }

    // swap x and xn
    temp = xn.vector;
    xn.vector = x.vector;
    x.vector = temp;

    iter++;
  } while(iter <= MAX_ITERATIONS);

  if (iter >= MAX_ITERATIONS)
    std::cout << "Max iteration reached" << std::endl;

  FreeVector(x);
  FreeVector(xn);
  FreeVector(c);
  FreeVector(r);
  FreeVector(e);
  FreeMatrix(m);
  FreeMatrix(dInv);
  FreeMatrix(lu);
}

void GaussSeidel(const Matrix& a, const Vector& b, const Vector& xStart) {
  Vector x;
  Vector xn;
  Vector c;
  Vector r;
  Vector e;
  Matrix m;
  Matrix lt;
  Matrix ltInv;
  Matrix u;

  int iter = 1;
  long double* temp;

  CreateVector(x, xStart.size);
  CreateVector(xn, xStart.size);
  CreateVector(c, xStart.size);
  CreateVector(r, xStart.size);
  CreateVector(e, xStart.size);
  CreateMatrix(m, a.size);
  CreateMatrix(lt, a.size);
  CreateMatrix(ltInv, a.size);
  CreateMatrix(u, a.size);

  // copy a to lt
  for (int i = 0; i < a.size; i++)
    for (int j = 0; j < a.size; j++)
      lt.matrix[i][j] = a.matrix[i][j];

  LUSplit(lt, u);
  LowerTriangularInverse(lt, ltInv);
  MMMultiplication(ltInv, u, m);
  MSMultiplication(m, -1.0L);

  MVMultiplication(ltInv, b, c);

  for (int i = 0; i < xStart.size; i++)
    x.vector[i] = xStart.vector[i];

  if (MatrixMaxNorm(m) >= 1.0L) {
    std::cout << "Method not converging" << std::endl;

    FreeVector(x);
    FreeVector(xn);
    FreeVector(c);
    FreeVector(r);
    FreeVector(e);
    FreeMatrix(m);
    FreeMatrix(lt);
    FreeMatrix(ltInv);
    FreeMatrix(u);

    return;
  }

  do {
    // calculate xn
    MVMultiplication(m, x, xn);
    VectorAdd(xn, c);

    // residuum
    MVMultiplication(a, xn, r);
    VectorSubtract(r, b);

    // err
    for (int i = 0; i < e.size; i++)
      e.vector[i] = std::fabs(xn.vector[i] - x.vector[i]);

    std::cout << "--- x" << iter << " ---" << std::endl;
    PrintVector(xn);

    std::cout << "--- residuum ---" << std::endl;
    PrintVector(r);
    std::cout << "MaxNorm: " << VectorMaxNorm(r) << std::endl;

    std::cout << "--- error estimator ---" << std::endl;
    PrintVector(e);
    std::cout << "MaxNorm: " << VectorMaxNorm(e) << std::endl;
    std::cout << std::endl;

    if (VectorMaxNorm(e) < TOLX && VectorMaxNorm(r) < TOLF) {
      std::cout << "Precision reached" << std::endl;
      break;
    }

    // swap x and xn
    temp = xn.vector;
    xn.vector = x.vector;
    x.vector = temp;

    iter++;
  } while(iter <= MAX_ITERATIONS);

  if (iter >= MAX_ITERATIONS)
    std::cout << "Max iteration reached" << std::endl;

  FreeVector(x);
  FreeVector(xn);
  FreeVector(c);
  FreeVector(r);
  FreeVector(e);
  FreeMatrix(m);
  FreeMatrix(lt);
  FreeMatrix(ltInv);
  FreeMatrix(u);
}

void SOR(const Matrix& a, const Vector& b, const Vector& xStart,
    long double omega) {
  Vector x;
  Vector xn;
  Vector c;
  Vector r;
  Vector e;
  Matrix m;
  Matrix lt;
  Matrix ltInv;
  Matrix u;

  int iter = 1;
  long double* temp;

  CreateVector(x, xStart.size);
  CreateVector(xn, xStart.size);
  CreateVector(c, xStart.size);
  CreateVector(r, xStart.size);
  CreateVector(e, xStart.size);
  CreateMatrix(m, a.size);
  CreateMatrix(lt, a.size);
  CreateMatrix(ltInv, a.size);
  CreateMatrix(u, a.size);

  // copy a to lt
  for (int i = 0; i < a.size; i++)
    for (int j = 0; j < a.size; j++)
      lt.matrix[i][j] = a.matrix[i][j];

  LUSplit(lt, u);

  for (int i = 0; i < lt.size; i++) {
    u.matrix[i][i] = (1 - 1 / omega) * lt.matrix[i][i];
    lt.matrix[i][i] = lt.matrix[i][i] / omega;
  }

  LowerTriangularInverse(lt, ltInv);
  MMMultiplication(ltInv, u, m);
  MSMultiplication(m, -1.0L);

  MVMultiplication(ltInv, b, c);

  for (int i = 0; i < xStart.size; i++)
    x.vector[i] = xStart.vector[i];

  if (MatrixMaxNorm(m) >= 1.0L || omega < 0 || omega > 2) {
    std::cout << "Method not converging" << std::endl;

    FreeVector(x);
    FreeVector(xn);
    FreeVector(c);
    FreeVector(r);
    FreeVector(e);
    FreeMatrix(m);
    FreeMatrix(lt);
    FreeMatrix(ltInv);
    FreeMatrix(u);

    return;
  }

  do {
    // calculate xn
    MVMultiplication(m, x, xn);
    VectorAdd(xn, c);

    // residuum
    MVMultiplication(a, xn, r);
    VectorSubtract(r, b);

    // err
    for (int i = 0; i < e.size; i++)
      e.vector[i] = std::fabs(xn.vector[i] - x.vector[i]);

    std::cout << "--- x" << iter << " ---" << std::endl;
    PrintVector(xn);

    std::cout << "--- residuum ---" << std::endl;
    PrintVector(r);
    std::cout << "MaxNorm: " << VectorMaxNorm(r) << std::endl;

    std::cout << "--- error estimator ---" << std::endl;
    PrintVector(e);
    std::cout << "MaxNorm: " << VectorMaxNorm(e) << std::endl;
    std::cout << std::endl;

    if (VectorMaxNorm(e) < TOLX && VectorMaxNorm(r) < TOLF) {
      std::cout << "Precision reached" << std::endl;
      break;
    }

    // swap x and xn
    temp = xn.vector;
    xn.vector = x.vector;
    x.vector = temp;

    iter++;
  } while(iter <= MAX_ITERATIONS);

  if (iter >= MAX_ITERATIONS)
    std::cout << "Max iteration reached" << std::endl;

  FreeVector(x);
  FreeVector(xn);
  FreeVector(c);
  FreeVector(r);
  FreeVector(e);
  FreeMatrix(m);
  FreeMatrix(lt);
  FreeMatrix(ltInv);
  FreeMatrix(u);
}
