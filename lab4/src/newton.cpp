#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "newton.hpp"

int main() {
  newton();
  return 0;
}

void newton() {
  int i = 0;

  std::vector<double> xs = { 1, 1, 1 };
  std::vector<double> xsn(3, 0.0);
  std::vector<double> fxs(3, 0.0);

  std::vector<double> es(3, 0.0);
  std::vector<double> residuums(3, 0.0);

  std::vector<std::vector<double>> jac;

  do {
    jac = jacobian(xs);
    jac = inverse(jac);

    fxs = f(xs);
    fxs = multiply(jac, fxs);

    for (int j = 0; j < 3; j++)
      xsn[j] = xs[j] - fxs[j];

    for (int j = 0; j < 3; j++)
      es[j] = xsn[j] - xs[j];

    residuums = f(xsn);

    std::cout << "xsn: ";
    print_vector(xsn);
    std::cout << "\t" << "es: ";
    print_vector(es);
    std::cout << "\t" << "res: ";
    print_vector(residuums);
    std::cout << std::endl;

    if (es[0] <= TOLX && es[1] <= TOLX && es[2] <= TOLX) {
      std::cout << "blad mniejszy od TOLX" << std::endl;
      return;
    }

    if (residuums[0] <= TOLF && residuums[1] <= TOLF 
        && residuums[2] <= TOLF) {
      std::cout << "residuum mniejsze od TOLF" << std::endl;
      return;
    }

    xs = xsn;
    i++;
  } while (i < MAX_ITERATIONS);
  
  std::cout << "maksymalna ilosc iteracji" << std::endl;
}

std::vector<double> f(std::vector<double> xs) {
  if (xs.size() != 3) throw std::runtime_error("invalid vector size");

  std::vector<double> result = {
    xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2] - 4,
    xs[0] * xs[0] + xs[1] * xs[1] * 0.5 - 1,
    xs[0] * xs[1] - 0.5
  };

  return result;
}

std::vector<std::vector<double>> jacobian(std::vector<double> xs) {
  if (xs.size() != 3) throw std::runtime_error("invalid vector size");

  std::vector<std::vector<double>> result = {
    { 2 * xs[0], 2 * xs[1], 2 * xs[2] },
    { 2 * xs[0], xs[1], 0},
    { xs[1], xs[0], 0 }
  };

  return result;
}

double det(std::vector<std::vector<double>> matrix) {
  if (matrix.size() != 3) throw std::runtime_error("invalid matrix size");

  for(auto& r : matrix)
    if (r.size() != 3)
      throw std::runtime_error("invalid matrix size");

  double m1 = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
  double m2 = matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0];
  double m3 = matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0];

  return matrix[0][0] * m1 - matrix[0][1] * m2 + matrix[0][2] * m3;
}

std::vector<std::vector<double>>
  adj(std::vector<std::vector<double>> matrix) {
  if (matrix.size() != 3) throw std::runtime_error("invalid matrix size");

  for(auto& r : matrix)
    if (r.size() != 3)
      throw std::runtime_error("invalid matrix size");

  std::vector<std::vector<double>> result(3, std::vector<double>(3));

  result[0][0] = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
  result[0][1] = matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2];
  result[0][2] = matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1];

  result[1][0] = matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2];
  result[1][1] = matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0];
  result[1][2] = matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2];

  result[2][0] = matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0];
  result[2][1] = matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1];
  result[2][2] = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  
  return result;
}

std::vector<std::vector<double>>
  inverse(std::vector<std::vector<double>> matrix) {
  if (matrix.size() != 3) throw std::runtime_error("invalid matrix size");

  for(auto& r : matrix)
    if (r.size() != 3)
      throw std::runtime_error("invalid matrix size");

  double d = det(matrix);

  if (d == 0) throw std::runtime_error("singular matrix");

  std::vector<std::vector<double>> a = adj(matrix);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      a[i][j] = a[i][j] / d;

  return a;
}

std::vector<double> multiply (std::vector<std::vector<double>> matrix,
  std::vector<double> xs) {
  if (xs.size() != 3) throw std::runtime_error("invalid vector size");
  if (matrix.size() != 3) throw std::runtime_error("invalid matrix size");

  for(auto& r : matrix)
    if (r.size() != 3)
      throw std::runtime_error("invalid matrix size");

  std::vector<double> result = { 0, 0, 0 };

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result[i] += matrix[i][j] * xs[j];

  return result;
}

void print_vector(std::vector<double> xs) {
  for (auto& x : xs)
    std::cout << x << ", ";
}
