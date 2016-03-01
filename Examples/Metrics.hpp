#ifndef CME212_METRICS_HPP
#define CME212_METRICS_HPP

/**
 * @file Metrics.hpp
 * @brief Template class to keep track of metric units
 */

#include "CME212/Point.hpp"

namespace Metrics {
/**
 * @class Unit
 * @brief Template class representing a metric unit
 *
 * The template parameters i, j, and k represent the exponents of the mass,
 * distance, and time of the unit, respectively.
 */
template <class T, int i, int j, int k>
struct Unit {
  using value_type = T;
  T value;

  // CONSTRUCTORS

  /* Unit with given value */
  constexpr Unit(const T& v = T())
      : value(v) {
  }

  /* Output stream operator */
  friend std::ostream& operator<<(std::ostream& s, const Unit& u) {
    return s << u.value << " [g^" << i << "m^" << j << "s^" << k << "]";
  }
};

// TYPEDEFS - Common units

using Scalar        = Unit<double, 0, 0,  0>;
using Mass          = Unit<double, 1, 0,  0>;
using Distance      = Unit<double, 0, 1,  0>;
using Distance3     = Unit<Point,  0, 1,  0>;
using Time          = Unit<double, 0, 0,  1>;
using Velocity      = Unit<double, 0, 1, -1>;
using Velocity3     = Unit<Point,  0, 1, -1>;
using Acceleration  = Unit<double, 0, 1, -2>;
using Acceleration3 = Unit<Point,  0, 1, -2>;
using Force         = Unit<double, 1, 1, -2>;
using Force3        = Unit<Point,  1, 1, -2>;
using SpringConst   = Unit<double, 1, 0, -2>;
using DampingConst  = Unit<double, 1, 0, -1>;

// USER DEFINED LITERALS -- Sugar for constructing units

constexpr Mass operator"" _g(long double v) {
  return {v};
}
constexpr Mass operator"" _g(unsigned long long v) {
  return {v};
}
constexpr Time operator"" _s(long double v) {
  return {v};
}
constexpr Time operator"" _s(unsigned long long v) {
  return {v};
}
constexpr Distance operator"" _m(long double v) {
  return {v};
}
constexpr Distance operator"" _m(unsigned long long v) {
  return {v};
}

// COMPARATORS

template <class T, int i, int j, int k>
inline
bool operator==(const Unit<T, i, j, k>& a,
                const Unit<T, i, j, k>& b) {
  return a.value == b.value;
}


// UNARY FUNCTIONS

/* Return @a -a */
template <class T, int i, int j, int k>
inline
Unit<T, i, j, k> operator-(const Unit<T, i, j, k>& a) {
  return {-a.value};
}

/* Return @a a */
template <class T, int i, int j, int k>
inline
const Unit<T, i, j, k>& operator+(const Unit<T, i, j, k>& a) {
  return a;
}

// ADDITION AND SUBTRACTION

/* Addition and subtraction of units should not change their type */
template <class T, class U, int i, int j, int k>
inline
Unit<T, i, j, k>& operator+=(Unit<T, i, j, k>& a,
                             const Unit<U, i, j, k>& b) {
  a.value += b.value;
  return a;
}

template <class T, class U, int i, int j, int k>
inline Unit<T, i, j, k>& operator-=(Unit<T, i, j, k>& a,
                                    const Unit<U, i, j, k>& b) {
  a.value -= b.value;
  return a;
}

// MULTIPLICATION AND DIVISION

/* Only scalars can scale a Unit */
template <class T, int i, int j, int k>
inline
Unit<T, i, j, k>& operator*=(Unit<T, i, j, k>& a,
                             const Scalar& b) {
  a.value *= b.value;
  return a;
}

/* Only scalars can scale a Unit */
template <class T, int i, int j, int k>
inline
Unit<T, i, j, k>& operator/=(Unit<T, i, j, k>& a,
                             const Scalar& b) {
  a.value /= b.value;
  return a;
}

/* Addition and subtraction units does not change the units of the result, but
 * may change the type of the result. The type of the value
 * is determined as the type of the sum and difference of those values.
 */

/** Define a helper class that is the sum of two types */
template <class T, class U>
using add_type = decltype(std::declval<T>() + std::declval<U>());

template <class T, class U, int i, int j, int k>
inline
Unit<add_type<T,U>, i, j, k> operator+(const Unit<T, i, j, k>& a,
                                       const Unit<U, i, j, k>& b) {
  return {a.value + b.value};
}

/** Define a helper class that is the difference of two types */
template <class T, class U>
using sub_type = decltype(std::declval<T>() - std::declval<U>());

template <class T, class U, int i, int j, int k>
inline
Unit<sub_type<T,U>, i, j, k> operator-(const Unit<T, i, j, k>& a,
                                       const Unit<U, i, j, k>& b) {
  return {a.value - b.value};
}

/* Multiplication and division of units changes the type of the result since
 * the exponents are added or subtracted, respectively.  The type of the value
 * is determined as the type of the products and divisions of those values.
 */

/** Define a helper type that is the product of two types */
template <class T, class U>
using prod_type = decltype(std::declval<T>() * std::declval<U>());

template <class T, class U, int i, int j, int k, int l, int m, int n>
inline Unit<prod_type<T,U>, i+l, j+m, k+n> operator*(const Unit<T, i, j, k>& a,
                                                     const Unit<U, l, m, n>& b) {
  return {a.value * b.value};
}

/** Define a helper type that is the division of two types */
template <class T, class U>
using div_type = decltype(std::declval<T>() / std::declval<U>());

template <class T, class U, int i, int j, int k, int l, int m, int n>
inline
Unit<div_type<T,U>, i-l, j-m, k-n> operator/(const Unit<T, i, j, k>& a,
                                             const Unit<U, l, m, n>& b) {
  return {a.value / b.value};
}

};

#endif // CME212_METRICS_HPP
