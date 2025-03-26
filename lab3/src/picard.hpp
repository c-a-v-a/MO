#ifndef PICARD_HPP
#define PICARD_HPP

static const int MAX_ITERATIONS = 50;
static const double TOLX = 0.000001;
static const double TOLF = 0.000001;

int main();

void picard(double (*f)(double), double (*phi)(double));

double phi1(double x);

double f1(double x);

double phi2(double x);

double f2(double x);

#endif
