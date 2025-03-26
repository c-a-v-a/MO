#ifndef SIECZNE_HPP
#define SIECZNE_HPP

static const int MAX_ITERATIONS = 50;
static const double TOLX = 0.000001;
static const double TOLF = 0.000001;

int main();

void sieczne(double (*f)(double), double x, double xn);

double f1(double x);

double f2(double x);

#endif
