#ifndef MAIN_H
#define MAIN_H

static const long double T_MAX = 1.0l;

static const long double first = 2.0l;

long double f(long double t, long double y);

long double real(long double t);

long double bme_step(long double t, long double y, long double dt);

long double bme(long double dt);

long double pme_step(long double t, long double y, long double dt);

long double pme(long double dt);

long double pmt_step(long double t, long double y, long double dt);

long double pmt(long double dt);

#endif
