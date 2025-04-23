#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>

#include "main.h"

int main() {
  double start = 1.0e-2;
  double end = 1.1e-40;
  double PI = std::numbers::pi;
  double h = start;
  std::ofstream doublefile("double.dat");

  do {
    doublefile 
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << std::log10(h) << " ";

    // x = 0
    doublefile
      << std::log10(std::fabs(df(0.0) - forward<double>(0.0, h)))
      << " " 
      << std::log10(std::fabs(df(0.0) - threepoint_left<double>(0.0, h)))
      << " ";

    // x = PI/4
    doublefile 
      << std::log10(std::fabs(df(PI / 4) - forward<double>(PI / 4, h)))
      << " "
      << std::log10(std::fabs(df(PI / 4) - backward<double>(PI / 4, h)))
      << " "
      << std::log10(std::fabs(df(PI / 4) - central<double>(PI / 4, h)))
      << " ";

    // x = PI/2
    doublefile
      << std::log10(std::fabs(df(PI / 2) - backward<double>(PI / 2, h)))
      << " " 
      << std::log10(std::fabs(df(PI / 2) - threepoint_right<double>(PI / 2, h)));

    h /= 10.0;

    doublefile << std::endl;
  } while (h > end);

  doublefile.close();

  long double lstart = 1.0e-2L;
  long double lend = 1.1e-40L;
  long double LPI = std::numbers::pi_v<long double>;
  long double lh = start;
  std::ofstream ldoublefile("ldouble.dat");

  do {
    ldoublefile 
      << std::setprecision(std::numeric_limits<long double>::max_digits10)
      << std::log10(lh) << " ";

    // x = 0
    ldoublefile
      << std::log10(std::fabs(df(0.0L) - forward<long double>(0.0L, lh)))
      << " " 
      << std::log10(std::fabs(df(0.0L) - threepoint_left<long double>(0.0L, lh)))
      << " ";

    // x = PI/4
    ldoublefile 
      << std::log10(std::fabs(df(LPI / 4) - forward<long double>(LPI / 4, lh)))
      << " "
      << std::log10(std::fabs(df(LPI / 4) - backward<long double>(LPI / 4, lh)))
      << " "
      << std::log10(std::fabs(df(LPI / 4) - central<long double>(LPI / 4, lh)))
      << " ";

    // x = PI/2
    ldoublefile
      << std::log10(std::fabs(df(LPI / 2) - backward<long double>(LPI / 2, lh)))
      << " " 
      << std::log10(std::fabs(df(LPI / 2) - threepoint_right<long double>(LPI / 2, lh)));

    lh /= 10.0L;

    ldoublefile << std::endl;
  } while (lh > lend);

  ldoublefile.close();
}
