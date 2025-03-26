#ifndef NEWTON_HPP
#define NEWTON_HPP

static const int MAX_ITERATIONS = 50;
static const double TOLX = 0.000001;
static const double TOLF = 0.000001;

int main();

void newton(double (*f)(double), double (*df)(double));

double f1(double x);

double df1(double x);

double f2(double x);

double df2(double x);

#endif
