#pragma once

#include <iosfwd>

#include "../rational.h"
#include "integer.h"

namespace poly {

  /** Wrapper for a rational. */
  class Rational {
    /** Actual rational number. */
    lp_rational_t mRat;

   public:
    /** Construct from an internal lp_rational_t pointer. */
    explicit Rational(const lp_rational_t* r);
    /** Construct as zero. */
    Rational();
    /** Construct from an int. */
    Rational(int i);
    /** Copy from a rational. */
    Rational(const Rational& r);
    /** Move from a rational. */
    Rational(Rational&& r);

    /** Copy from numerator and denominator. */
    Rational(long num, unsigned long denom);
    /** Copy from numerator and denominator. */
    Rational(const Integer& num, const Integer& denom);
    /** Copy from Integer. */
    Rational(const Integer& i);
    /** Construct as d. */
    Rational(double d);

    /** Custom destructor. */
    ~Rational();

    /** Copy from a rational. */
    Rational& operator=(const Rational& r);
    /** Move from a rational. */
    Rational& operator=(Rational&& r);

    /** Get a non-const pointer to the internal lp_rational_t. Handle with care!
     */
    lp_rational_t* get_internal();
    /** Get a const pointer to the internal lp_rational_t. */
    const lp_rational_t* get_internal() const;
  };

  /** Make sure that we can cast between Rational and lp_rational_t. */
  static_assert(sizeof(Rational) == sizeof(lp_rational_t),
                "Please check the size of Rational.");
  namespace detail {
    /** Non-const cast from a Rational to a lp_rational_t. */
    inline lp_rational_t* cast_to(Rational* i) {
      return reinterpret_cast<lp_rational_t*>(i);
    }
    /** Const cast from a Rational to a lp_rational_t. */
    inline const lp_rational_t* cast_to(const Rational* i) {
      return reinterpret_cast<const lp_rational_t*>(i);
    }
    /** Non-const cast from a lp_rational_t to a Rational. */
    inline Rational* cast_from(lp_rational_t* i) {
      return reinterpret_cast<Rational*>(i);
    }
    /** Const cast from a lp_rational_t to a Rational. */
    inline const Rational* cast_from(const lp_rational_t* i) {
      return reinterpret_cast<const Rational*>(i);
    }
  }  // namespace detail

  /** Stream the given Rational to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Rational& r);

  /** Give a double approximation. */
  double to_double(const Rational& r);

  /** Give the sign of a rational. */
  int sgn(const Rational& r);

  /** Compare two rational. */
  bool operator==(const Rational& lhs, const Rational& rhs);
  /** Compare two rational. */
  bool operator!=(const Rational& lhs, const Rational& rhs);
  /** Compare two rational. */
  bool operator<(const Rational& lhs, const Rational& rhs);
  /** Compare two rational. */
  bool operator<=(const Rational& lhs, const Rational& rhs);
  /** Compare two rational. */
  bool operator>(const Rational& lhs, const Rational& rhs);
  /** Compare two rational. */
  bool operator>=(const Rational& lhs, const Rational& rhs);

  /** Compare a rational and an integer. */
  bool operator==(const Rational& lhs, const Integer& rhs);
  /** Compare a rational and an integer. */
  bool operator!=(const Rational& lhs, const Integer& rhs);
  /** Compare a rational and an integer. */
  bool operator<(const Rational& lhs, const Integer& rhs);
  /** Compare a rational and an integer. */
  bool operator<=(const Rational& lhs, const Integer& rhs);
  /** Compare a rational and an integer. */
  bool operator>(const Rational& lhs, const Integer& rhs);
  /** Compare a rational and an integer. */
  bool operator>=(const Rational& lhs, const Integer& rhs);

  /** Compare an integer and a rational. */
  bool operator==(const Integer& lhs, const Rational& rhs);
  /** Compare an integer and a rational. */
  bool operator!=(const Integer& lhs, const Rational& rhs);
  /** Compare an integer and a rational. */
  bool operator<(const Integer& lhs, const Rational& rhs);
  /** Compare an integer and a rational. */
  bool operator<=(const Integer& lhs, const Rational& rhs);
  /** Compare an integer and a rational. */
  bool operator>(const Integer& lhs, const Rational& rhs);
  /** Compare an integer and a rational. */
  bool operator>=(const Integer& lhs, const Rational& rhs);

  /** Swap two rationals. */
  void swap(Rational& lhs, Rational& rhs);

  /** Add two rationals. */
  Rational operator+(const Rational& lhs, const Rational& rhs);
  /** Add a rational and an integer. */
  Rational operator+(const Rational& lhs, const Integer& rhs);
  /** Add an integer and a rational. */
  Rational operator+(const Integer& lhs, const Rational& rhs);

  /** Subtract two rationals. */
  Rational operator-(const Rational& lhs, const Rational& rhs);
  /** Negate a rational. */
  Rational operator-(const Rational& r);

  /** Invert a rational. */
  Rational inverse(const Rational& r);

  /** Multiply two rationals. */
  Rational operator*(const Rational& lhs, const Rational& rhs);
  /** Compute lhs * 2^n. */
  Rational mul_2exp(const Rational& lhs, unsigned n);
  /** Compute r^n. */
  Rational pow(const Rational& r, unsigned n);

  /** Divide two rationals. */
  Rational operator/(const Rational& lhs, const Rational& rhs);
  /** Compute lhs / 2^n. */
  Rational div_2exp(const Rational& lhs, unsigned n);

  /** Return the numerator of a rational. */
  const Integer& numerator(const Rational& r);
  /** Return the denominator of a rational. */
  const Integer& denominator(const Rational& r);

  /** Check if a rational is integral. */
  bool is_integer(const Rational& r);

  /** Compute the ceiling of a rational. */
  Integer ceil(const Rational& r);
  /** Compute the floor of a rational. */
  Integer floor(const Rational& r);

}  // namespace poly
