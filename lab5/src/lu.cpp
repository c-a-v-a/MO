#include <iostream>
#include <stdexcept>
#include <vector>

#include "lu.hpp"

bool IsSquare(const Matrix& x) {
  int size = x.size();

  for (auto y : x)
    if (y.size() != size) return false;

  return true;
}

Vector GetColumn(const Matrix& x, int column) {
  Vector v;

  for (auto y : x) {
    if (y.size() < column)
      throw std::runtime_error("Invalid column index");

    v.push_back(y.at(column));
  }

  return v;
}

int MaxValueIndex(const Vector& x, int start) {
  int maxIndex = 0;

  for (int i = start; i < x.size(); i++)
    if (x.at(maxIndex) < x.at(i))
      maxIndex = i;

  return maxIndex;
}

void LUFactorization(Matrix& l, Matrix& u, Vector& p) {
  if (!IsSquare(u))
    throw std::runtime_error("Matrix not square");

  ZeroMatrix(l, u.size());
  PermutationVector(p, u.size());

  for (int i = 0; i < u.size(); i++) {
    int maxIndex = MaxValueIndex(GetColumn(u, i), i+1);

    if (maxIndex > i)
      std::swap(p.at(i), p.at(maxIndex));

    int topRow = p.at(i);

    l.at(topRow).at(i) = 1.0L;

    if (u.at(topRow).at(i) == 0)
      topRow = 0;
    
    for (int j = i + 1; j < u.size(); j++) {
      int row = p.at(j);
      long double k = u.at(row).at(i) / u.at(topRow).at(i);

      l.at(row).at(i) = k;
      VectorSubtract(u.at(row), u.at(topRow), k);
    }
  }
}

void PrintVector(const Vector& x) {
  for (auto y : x)
    std::cout << y << " ";

  std::cout << std::endl;
}

void PrintMatrix(const Matrix& x) {
  for (auto y : x)
    PrintVector(y);
}

void PrintPVector(const Vector& x, const Vector& p) {
  for (int i = 0; i < x.size(); i++)
    std::cout << x.at(p.at(i)) << " ";

  std::cout << std::endl;
}

void PrintPMatrix(const Matrix& x, const Vector& p) {
  for (int i = 0; i < x.size(); i++) {
    Vector y = x.at(p.at(i));

    for (int j = 0; j < y.size(); j++)
      std::cout << y.at(j) << " ";

    std::cout << std::endl;
  }
}

void ZeroMatrix(Matrix& x, int size) {
  x.resize(0);

  for (int i = 0; i < size; i++) {
    Vector v(size, 0.0L);

    x.push_back(v);
  }
}

void PermutationVector(Vector& p, int size) {
  p.resize(0);

  for (int i = 0; i < size; i++)
    p.push_back(i);
}

void VectorSubtract(Vector& x, const Vector& y, long double k /*= 1.0L*/) {
  if (x.size() != y.size())
    throw std::runtime_error("Vectors not the same size");

  for (int i = 0; i < x.size(); i++) {
    x.at(i) -= k * y.at(i);
  }
}

void CalculateY(const Matrix& l, Vector& y, const Vector& b, const Vector& p) {
  if (!IsSquare(l))
    throw std::runtime_error("Matrix L is not square");
  if (l.size() != b.size())
    throw std::runtime_error("Vector b is wrong size");

  Vector x(l.size(), 0.0L);

  x.at(p.at(0)) = b.at(p.at(0)) / l.at(p.at(0)).at(0);

  for (int i = 1; i < l.size(); i++) {
    long double sum = b.at(p.at(i));

    for (int j = 0; j < i; j++) {
      sum -= l.at(p.at(i)).at(j) * x.at(p.at(j));
    }

    x.at(p.at(i)) = sum / l.at(p.at(i)).at(i);
  }

  y = x;
}

void CalculateX(const Matrix& u, Vector& x, const Vector& y, const Vector& p) {
   if (!IsSquare(u))
    throw std::runtime_error("Matrix U is not square");
  if (u.size() != y.size())
    throw std::runtime_error("Vector y is wrong size"); 

  Vector z(u.size(), 0.0L);

  z.at(p.back()) = y.at(p.back()) / u.at(p.back()).back();

  for (int i = u.size() - 2; i >= 0; i--) {
    long double sum = y.at(p.at(i));

    for (int j = i + 1; j < u.size(); j++) {
      sum -= u.at(p.at(i)).at(j) * z.at(p.at(j));
    }

    z.at(p.at(i)) = sum / u.at(p.at(i)).at(i);
  }

  x = z;
}

void Solve(Matrix& a, Vector& x, const Vector& b) {
  if (!IsSquare(a))
    throw std::runtime_error("Matrix A is not square");
  if(a.size() != b.size())
    throw std::runtime_error("Vector b is wrong size"); 

  #ifndef NOPRINT
  std::cout << "=== SOLVE ===" << std::endl;
  std::cout << "--- MATRIX A ---" << std::endl;
  PrintMatrix(a);
  std::cout << std::endl;
  std::cout << "--- VECTOR B ---" << std::endl;
  PrintVector(b);
  std::cout << std::endl;
  #endif

  Vector p;
  Vector y;
  Matrix l;

  LUFactorization(l, a, p);

  #ifndef NOPRINT
  std::cout << "--- MATRIX L ---" << std::endl;
  PrintPMatrix(l, p);
  std::cout << std::endl;
  std::cout << "--- MATRIX U ---" << std::endl;
  PrintPMatrix(a, p);
  std::cout << std::endl;
  std::cout << "--- PERMUTATION VECTOR ---" << std::endl;
  PrintVector(p);
  std::cout << std::endl;
  #endif

  CalculateY(l, y, b, p);

  #ifndef NOPRINT
  std::cout << "--- VECTOR Y ---" << std::endl;
  PrintPVector(y, p);
  std::cout << std::endl;
  #endif

  CalculateX(a, x, y, p);

  #ifndef NOPRINT
  std::cout << "--- VECTOR X ---" << std::endl;
  PrintPVector(x, p);
  std::cout << std::endl;
  #endif

  Vector z(x.size(), 0.0L);

  for (int i = 0; i < z.size(); i++)
    z.at(i) = x.at(p.at(i));

  x = z;
}
