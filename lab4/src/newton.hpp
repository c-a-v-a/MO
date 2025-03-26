#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <vector>

static const int MAX_ITERATIONS = 50;
static const double TOLX = 0.000000001;
static const double TOLF = 0.000000001;

int main();

void newton();

std::vector<double> f(std::vector<double> xs);

std::vector<std::vector<double>> jacobian(std::vector<double> xs);

double det(std::vector<std::vector<double>> matrix);

std::vector<std::vector<double>>
  adj(std::vector<std::vector<double>> matrix);

std::vector<std::vector<double>>
  inverse(std::vector<std::vector<double>> matrix);

std::vector<double> multiply (std::vector<std::vector<double>> matrix,
  std::vector<double> xs);

void print_vector(std::vector<double> xs);

#endif
