#include <cmath>
#include <fstream>
#include <iostream>

#include "main.h"

int main() {
  long double dt = 0.1l;
  long double error = 0.0l;

  error = bme(dt);
  std::cout << error << std::endl;

  error = pme(dt);
  std::cout << error << std::endl;

  error = pmt(dt);
  std::cout << error << std::endl;
}

long double f(long double t, long double y) {
  return ((100 * t + 10) / (t + 1)) * (y - 1);
}

long double real(long double t) {
  return 1 + std::powl((1 + t), 90.0l) * std::exp(-100.0l * t);
}

long double bme_step(long double t, long double y, long double dt) {
  return y + f(t, y) * dt;
}

long double bme(long double dt) {
  long double t = 0.0l;
  long double y = first;
  long double error = 0.0;
  
  while (t <= T_MAX) {
    std::cout << "T: " << t << ", Y: " << y << std::endl;

    y = bme_step(t, y, dt);

    if (std::fabs(y - real(t)) > error)
      error = std::fabs(y - real(t));

    t += dt;
  }

  return error;
}

long double pme_step(long double t, long double y, long double dt) {
  long double newy = (1 / dt * y + (100.0l * (t + dt) + 10.0l) / (t + dt + 1.0l));
  newy /= (1 / dt + (100.0l * (t + dt) + 10.0l) / (t + dt + 1.0l));

  return newy;
}

long double pme(long double dt) {
  long double t = 0.0l;
  long double y = first;
  long double error = 0.0;
  
  while (t <= T_MAX) {
    std::cout << "T: " << t << ", Y: " << y << std::endl;

    y = pme_step(t, y, dt);

    if (std::fabs(y - real(t)) > error)
      error = std::fabs(y - real(t));

    t += dt;
  }

  return error;
}

long double pmt_step(long double t, long double y, long double dt) {
  long double newy = 1 / dt * y;
  newy += 0.5l * f(t, y) - (50.0l * (t + dt) + 5.0l) / (t + dt + 1.0l);
  newy /= (1 / dt - (50.0l * (t + dt) + 5.0l) / (t + dt + 1.0l));

  return newy;
}

long double pmt(long double dt) {
  long double t = 0.0l;
  long double y = first;
  long double error = 0.0;
  
  while (t <= T_MAX) {
    std::cout << "T: " << t << ", Y: " << y << std::endl;

    y = pmt_step(t, y, dt);

    if (std::fabs(y - real(t)) > error)
      error = std::fabs(y - real(t));

    t += dt;
  }

  return error;
}
