#pragma once

#include <iosfwd>

#include "../dyadic_rational.h"
#include "integer.h"
#include "rational.h"

namespace poly {

  /**
   * Implements a dyadic rational, wrapping lp_dyadic_rational_t.
   */
  class DyadicRational {
    /** The actual number. */
    lp_dyadic_rational_t mDRat;

   public:
    /** Construct as zero. */
    DyadicRational();
    /** Copy from a DyadicRational. */
    DyadicRational(const DyadicRational& dr);
    /** Move from a DyadicRational. */
    DyadicRational(DyadicRational&& dr);

    /** Construct as a / 2^n. */
    DyadicRational(long a, unsigned long n);
    /** Construct from an Integer. */
    explicit DyadicRational(const Integer& i);
    /** Construct from a double. */
    explicit DyadicRational(double d);
    /** Construct from an int. */
    DyadicRational(int i);
    /** Construct from a long. */
    DyadicRational(long i);

    /** Construct from an internal lp_dyadic_rational_t pointer. */
    explicit DyadicRational(const lp_dyadic_rational_t* dr);

    /** Custom destructor. */
    ~DyadicRational();

    /** Copy from a dr. */
    DyadicRational& operator=(const DyadicRational& dr);
    /** Move from a DyadicRational. */
    DyadicRational& operator=(DyadicRational&& dr);

    /** Cast to a regular Rational. */
    operator Rational() const;

    /** Get a non-const pointer to the internal lp_dyadic_rational_t. Handle
     * with care! */
    lp_dyadic_rational_t* get_internal();
    /** Get a const pointer to the internal lp_dyadic_rational_t. */
    const lp_dyadic_rational_t* get_internal() const;
  };

  /** Make sure that we can cast between DyadicRational and
   * lp_dyadic_rational_t. */
  static_assert(sizeof(DyadicRational) == sizeof(lp_dyadic_rational_t),
                "Please check the size of DyadicRational.");
  namespace detail {
    /** Non-const cast from an DyadicRational to a lp_dyadic_rational_t. */
    inline lp_dyadic_rational_t* cast_to(DyadicRational* i) {
      return reinterpret_cast<lp_dyadic_rational_t*>(i);
    }
    /** Const cast from an DyadicRational to a lp_dyadic_rational_t. */
    inline const lp_dyadic_rational_t* cast_to(const DyadicRational* i) {
      return reinterpret_cast<const lp_dyadic_rational_t*>(i);
    }
    /** Non-const cast from a lp_dyadic_rational_t to an DyadicRational. */
    inline DyadicRational* cast_from(lp_dyadic_rational_t* i) {
      return reinterpret_cast<DyadicRational*>(i);
    }
    /** Const cast from a lp_dyadic_rational_t to an DyadicRational. */
    inline const DyadicRational* cast_from(const lp_dyadic_rational_t* i) {
      return reinterpret_cast<const DyadicRational*>(i);
    }
  }  // namespace detail

  /** Stream the given DyadicRational to an output stream. */
  std::ostream& operator<<(std::ostream& os, const DyadicRational& dr);

  /** Get a double representation of a DyadicRational. */
  double to_double(const DyadicRational& dr);

  /** Get the sign of a DyadicRational. */
  int sgn(const DyadicRational& dr);

  /** Compare two DyadicRationals. */
  bool operator==(const DyadicRational& lhs, const DyadicRational& rhs);
  /** Compare two DyadicRationals. */
  bool operator!=(const DyadicRational& lhs, const DyadicRational& rhs);
  /** Compare two DyadicRationals. */
  bool operator<(const DyadicRational& lhs, const DyadicRational& rhs);
  /** Compare two DyadicRationals. */
  bool operator<=(const DyadicRational& lhs, const DyadicRational& rhs);
  /** Compare two DyadicRationals. */
  bool operator>(const DyadicRational& lhs, const DyadicRational& rhs);
  /** Compare two DyadicRationals. */
  bool operator>=(const DyadicRational& lhs, const DyadicRational& rhs);

  /** Compare a DyadicRational with an Integer. */
  bool operator==(const DyadicRational& lhs, const Integer& rhs);
  /** Compare a DyadicRational with an Integer. */
  bool operator!=(const DyadicRational& lhs, const Integer& rhs);
  /** Compare a DyadicRational with an Integer. */
  bool operator<(const DyadicRational& lhs, const Integer& rhs);
  /** Compare a DyadicRational with an Integer. */
  bool operator<=(const DyadicRational& lhs, const Integer& rhs);
  /** Compare a DyadicRational with an Integer. */
  bool operator>(const DyadicRational& lhs, const Integer& rhs);
  /** Compare a DyadicRational with an Integer. */
  bool operator>=(const DyadicRational& lhs, const Integer& rhs);

  /** Compare an Integer with a DyadicRational. */
  bool operator==(const Integer& lhs, const DyadicRational& rhs);
  /** Compare an Integer with a DyadicRational. */
  bool operator!=(const Integer& lhs, const DyadicRational& rhs);
  /** Compare an Integer with a DyadicRational. */
  bool operator<(const Integer& lhs, const DyadicRational& rhs);
  /** Compare an Integer with a DyadicRational. */
  bool operator<=(const Integer& lhs, const DyadicRational& rhs);
  /** Compare an Integer with a DyadicRational. */
  bool operator>(const Integer& lhs, const DyadicRational& rhs);
  /** Compare an Integer with a DyadicRational. */
  bool operator>=(const Integer& lhs, const DyadicRational& rhs);

  /** Compare a DyadicRational with a Rational. */
  bool operator==(const DyadicRational& lhs, const Rational& rhs);
  /** Compare a DyadicRational with a Rational. */
  bool operator!=(const DyadicRational& lhs, const Rational& rhs);
  /** Compare a DyadicRational with a Rational. */
  bool operator<(const DyadicRational& lhs, const Rational& rhs);
  /** Compare a DyadicRational with a Rational. */
  bool operator<=(const DyadicRational& lhs, const Rational& rhs);
  /** Compare a DyadicRational with a Rational. */
  bool operator>(const DyadicRational& lhs, const Rational& rhs);
  /** Compare a DyadicRational with a Rational. */
  bool operator>=(const DyadicRational& lhs, const Rational& rhs);

  /** Compare a Rational with a DyadicRational. */
  bool operator==(const Rational& lhs, const DyadicRational& rhs);
  /** Compare a Rational with a DyadicRational. */
  bool operator!=(const Rational& lhs, const DyadicRational& rhs);
  /** Compare a Rational with a DyadicRational. */
  bool operator<(const Rational& lhs, const DyadicRational& rhs);
  /** Compare a Rational with a DyadicRational. */
  bool operator<=(const Rational& lhs, const DyadicRational& rhs);
  /** Compare a Rational with a DyadicRational. */
  bool operator>(const Rational& lhs, const DyadicRational& rhs);
  /** Compare a Rational with a DyadicRational. */
  bool operator>=(const Rational& lhs, const DyadicRational& rhs);

  /** Swap two DyadicRationals. */
  void swap(DyadicRational& lhs, DyadicRational& rhs);

  /** Add two DyadicRationals. */
  DyadicRational operator+(const DyadicRational& lhs,
                           const DyadicRational& rhs);
  /** Add a DyadicRational and an Integer. */
  DyadicRational operator+(const DyadicRational& lhs, const Integer& rhs);
  /** Add an Integer and a DyadicRational. */
  DyadicRational operator+(const Integer& lhs, const DyadicRational& rhs);

  /** Subtract two DyadicRationals. */
  DyadicRational operator-(const DyadicRational& lhs,
                           const DyadicRational& rhs);

  /** Negate a DyadicRational. */
  DyadicRational operator-(const DyadicRational& dr);

  /** Multiply two DyadicRationals. */
  DyadicRational operator*(const DyadicRational& lhs,
                           const DyadicRational& rhs);

  /** Compute lhs * 2^n. */
  DyadicRational mul_2exp(const DyadicRational& lhs, unsigned long n);
  /** Compute lhs^n. */
  DyadicRational pow(const DyadicRational& lhs, unsigned long n);
  /** Compute lhs / 2^n. */
  DyadicRational div_2exp(const DyadicRational& lhs, unsigned long n);

  /** Get the numerator of a DyadicRational. */
  Integer numerator(const DyadicRational& r);
  /** Get the denominator of a DyadicRational. */
  Integer denominator(const DyadicRational& r);

  /** Check if a DyadicRational is an integer. */
  bool is_integer(const DyadicRational& r);

  /** Compute the ceiling of a DyadicRational. */
  Integer ceil(const DyadicRational& r);
  /** Compute the floor of a DyadicRational. */
  Integer floor(const DyadicRational& r);
}  // namespace poly
