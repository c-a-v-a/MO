#ifndef BISEKCJA_HPP
#define BISEKCJA_HPP

static const int MAX_ITERATIONS = 50;
static const double TOLX = 0.000001;
static const double TOLF = 0.000001;

int main();

void bisekcja(double (*f)(double), double a, double b);

double f1(double x);

double f2(double x);

#endif
