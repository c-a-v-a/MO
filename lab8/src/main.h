#ifndef MAIN_H
#define MAIN_H

#include <cmath>

int main();

template<typename T>
T f(T x) {
  return std::cos(x);
}

template<typename T>
T df(T x) {
  return (T)-1.0 * std::sin(x);
}

template<typename T>
T forward(T x, T h) {
  return (f(x + h) - f(x)) / h;
}

template<typename T>
T backward(T x, T h) {
  return (f(x) - f(x - h)) / h;
}

template<typename T>
T central(T x, T h) {
  return (f(x + h) - f(x - h)) / (2 * h);
}

template<typename T>
T threepoint_left(T x, T h) {
  return ((T)-1.5 * f(x) + (T)2.0 * f(x + h) - (T)0.5 * f(x + 2 * h)) / h;
}

template<typename T>
T threepoint_right(T x, T h) {
  return ((T)0.5 * f(x - 2 * h) - (T)2 * f(x - h) + (T)1.5 * f(x)) / h;
}

#endif
