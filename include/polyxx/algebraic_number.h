#pragma once

#include <iosfwd>

#include "../algebraic_number.h"
#include "dyadic_interval.h"
#include "dyadic_rational.h"
#include "integer.h"
#include "upolynomial.h"

namespace poly {

  /**
   * Implements a wrapper for lp_algebraic_number_t.
   */
  class AlgebraicNumber {
    /** The actual algebraic number. */
    lp_algebraic_number_t mValue;

   public:
    /** Copy from an internal lp_algebraic_number_t pointer. */
    explicit AlgebraicNumber(const lp_algebraic_number_t* an);
    /** Construct as zero. */
    AlgebraicNumber();
    /** Copy from the given AlgebraicNumber. */
    AlgebraicNumber(const AlgebraicNumber& an);
    /** Move from the given AlgebraicNumber. */
    AlgebraicNumber(AlgebraicNumber&& an);

    /** Construct from a DyadicRational */
    AlgebraicNumber(const DyadicRational& dr);

    /** Construct from a defining polynomial and an isolating interval. */
    AlgebraicNumber(UPolynomial&& poly, const DyadicInterval& di);
    /** Construct from a defining polynomial and an isolating interval. */
    AlgebraicNumber(const UPolynomial& poly, const DyadicInterval& di);
    /** Custom destructor. */
    ~AlgebraicNumber();
    /** Copy from the given AlgebraicNumber. */
    AlgebraicNumber& operator=(const AlgebraicNumber& an);
    /** Move from the given AlgebraicNumber. */
    AlgebraicNumber& operator=(AlgebraicNumber&& an);

    /** Get a non-const pointer to the internal lp_algebraic_number_t. Handle
     * with care! */
    lp_algebraic_number_t* get_internal();
    /** Get a const pointer to the internal lp_algebraic_number_t. */
    const lp_algebraic_number_t* get_internal() const;
  };

  /** Make sure that we can cast between AlgebraicNumber and
   * lp_algebraic_number_t. */
  static_assert(sizeof(AlgebraicNumber) == sizeof(lp_algebraic_number_t),
                "Please check the size of AlgebraicNumber.");
  namespace detail {
    /** Non-const cast from an AlgebraicNumber to a lp_algebraic_number_t. */
    inline lp_algebraic_number_t* cast_to(AlgebraicNumber* i) {
      return reinterpret_cast<lp_algebraic_number_t*>(i);
    }
    /** Const cast from an AlgebraicNumber to a lp_algebraic_number_t. */
    inline const lp_algebraic_number_t* cast_to(const AlgebraicNumber* i) {
      return reinterpret_cast<const lp_algebraic_number_t*>(i);
    }
    /** Non-const cast from a lp_algebraic_number_t to an AlgebraicNumber. */
    inline AlgebraicNumber* cast_from(lp_algebraic_number_t* i) {
      return reinterpret_cast<AlgebraicNumber*>(i);
    }
    /** Const cast from a lp_algebraic_number_t to an AlgebraicNumber. */
    inline const AlgebraicNumber* cast_from(const lp_algebraic_number_t* i) {
      return reinterpret_cast<const AlgebraicNumber*>(i);
    }
  }  // namespace detail

  /** Stream the given AlgebraicNumber to an output stream. */
  std::ostream& operator<<(std::ostream& os, const AlgebraicNumber& an);

  /** Swap two algebraic numbers. */
  void swap(AlgebraicNumber& lhs, AlgebraicNumber& rhs);
  /** Compute the sign of an algebraic number. */
  int sgn(const AlgebraicNumber& an);

  /** Compare two algebraic numbers. */
  bool operator==(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);
  /** Compare two algebraic numbers. */
  bool operator!=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);
  /** Compare two algebraic numbers. */
  bool operator<(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);
  /** Compare two algebraic numbers. */
  bool operator<=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);
  /** Compare two algebraic numbers. */
  bool operator>(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);
  /** Compare two algebraic numbers. */
  bool operator>=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs);

  /** Compare an algebraic number and an integer. */
  bool operator==(const AlgebraicNumber& lhs, const Integer& rhs);
  /** Compare an algebraic number and an integer. */
  bool operator!=(const AlgebraicNumber& lhs, const Integer& rhs);
  /** Compare an algebraic number and an integer. */
  bool operator<(const AlgebraicNumber& lhs, const Integer& rhs);
  /** Compare an algebraic number and an integer. */
  bool operator<=(const AlgebraicNumber& lhs, const Integer& rhs);
  /** Compare an algebraic number and an integer. */
  bool operator>(const AlgebraicNumber& lhs, const Integer& rhs);
  /** Compare an algebraic number and an integer. */
  bool operator>=(const AlgebraicNumber& lhs, const Integer& rhs);

  /** Compare an integer and an algebraic number. */
  bool operator==(const Integer& lhs, const AlgebraicNumber& rhs);
  /** Compare an integer and an algebraic number. */
  bool operator!=(const Integer& lhs, const AlgebraicNumber& rhs);
  /** Compare an integer and an algebraic number. */
  bool operator<(const Integer& lhs, const AlgebraicNumber& rhs);
  /** Compare an integer and an algebraic number. */
  bool operator<=(const Integer& lhs, const AlgebraicNumber& rhs);
  /** Compare an integer and an algebraic number. */
  bool operator>(const Integer& lhs, const AlgebraicNumber& rhs);
  /** Compare an integer and an algebraic number. */
  bool operator>=(const Integer& lhs, const AlgebraicNumber& rhs);

  /** Compare an algebraic number and a dyadic rational. */
  bool operator==(const AlgebraicNumber& lhs, const DyadicRational& rhs);
  /** Compare an algebraic number and a dyadic rational. */
  bool operator!=(const AlgebraicNumber& lhs, const DyadicRational& rhs);
  /** Compare an algebraic number and a dyadic rational. */
  bool operator<(const AlgebraicNumber& lhs, const DyadicRational& rhs);
  /** Compare an algebraic number and a dyadic rational. */
  bool operator<=(const AlgebraicNumber& lhs, const DyadicRational& rhs);
  /** Compare an algebraic number and a dyadic rational. */
  bool operator>(const AlgebraicNumber& lhs, const DyadicRational& rhs);
  /** Compare an algebraic number and a dyadic rational. */
  bool operator>=(const AlgebraicNumber& lhs, const DyadicRational& rhs);

  /** Compare a dyadic rational and an algebraic number. */
  bool operator==(const DyadicRational& lhs, const AlgebraicNumber& rhs);
  /** Compare a dyadic rational and an algebraic number. */
  bool operator!=(const DyadicRational& lhs, const AlgebraicNumber& rhs);
  /** Compare a dyadic rational and an algebraic number. */
  bool operator<(const DyadicRational& lhs, const AlgebraicNumber& rhs);
  /** Compare a dyadic rational and an algebraic number. */
  bool operator<=(const DyadicRational& lhs, const AlgebraicNumber& rhs);
  /** Compare a dyadic rational and an algebraic number. */
  bool operator>(const DyadicRational& lhs, const AlgebraicNumber& rhs);
  /** Compare a dyadic rational and an algebraic number. */
  bool operator>=(const DyadicRational& lhs, const AlgebraicNumber& rhs);

  /** Compare an algebraic number and a rational. */
  bool operator==(const AlgebraicNumber& lhs, const Rational& rhs);
  /** Compare an algebraic number and a rational. */
  bool operator!=(const AlgebraicNumber& lhs, const Rational& rhs);
  /** Compare an algebraic number and a rational. */
  bool operator<(const AlgebraicNumber& lhs, const Rational& rhs);
  /** Compare an algebraic number and a rational. */
  bool operator<=(const AlgebraicNumber& lhs, const Rational& rhs);
  /** Compare an algebraic number and a rational. */
  bool operator>(const AlgebraicNumber& lhs, const Rational& rhs);
  /** Compare an algebraic number and a rational. */
  bool operator>=(const AlgebraicNumber& lhs, const Rational& rhs);

  /** Compare a rational and an algebraic number. */
  bool operator==(const Rational& lhs, const AlgebraicNumber& rhs);
  /** Compare a rational and an algebraic number. */
  bool operator!=(const Rational& lhs, const AlgebraicNumber& rhs);
  /** Compare a rational and an algebraic number. */
  bool operator<(const Rational& lhs, const AlgebraicNumber& rhs);
  /** Compare a rational and an algebraic number. */
  bool operator<=(const Rational& lhs, const AlgebraicNumber& rhs);
  /** Compare a rational and an algebraic number. */
  bool operator>(const Rational& lhs, const AlgebraicNumber& rhs);
  /** Compare a rational and an algebraic number. */
  bool operator>=(const Rational& lhs, const AlgebraicNumber& rhs);

  /** Give a double approximation. */
  double to_double(const AlgebraicNumber& an);
  /** Give a rational approximation. */
  Rational to_rational_approximation(const AlgebraicNumber& an);

  /** Get lower bound of the isolating interval. */
  const DyadicRational& get_lower_bound(const AlgebraicNumber& an);
  /** Get upper bound of the isolating interval. */
  const DyadicRational& get_upper_bound(const AlgebraicNumber& an);
  /** Get the midpoint of the isolating interval. */
  DyadicRational midpoint_dyadic(const AlgebraicNumber& an);
  /** Get the midpoint of the isolating interval. */
  Rational midpoint_rational(const AlgebraicNumber& an);
  /** Refine the isolating interval. */
  void refine(AlgebraicNumber& an);
  /** Refine the isolating interval of a const algebraic number. */
  void refine_const(const AlgebraicNumber& an);

  /** Add two algebraic numbers. */
  AlgebraicNumber operator+(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs);
  /** Subtract two algebraic numbers. */
  AlgebraicNumber operator-(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs);
  /** Negate a algebraic number. */
  AlgebraicNumber operator-(const AlgebraicNumber& an);
  /** Multiply two algebraic numbers. */
  AlgebraicNumber operator*(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs);
  /** Compute lhs^n. */
  AlgebraicNumber pow(const AlgebraicNumber& lhs, unsigned n);

  /** Check whether an algebraic number is a rational.
   * This check is not complete.
   */
  bool is_rational(const AlgebraicNumber& an);
  /** Check whether an algebraic number is an integer. */
  bool is_integer(const AlgebraicNumber& an);
  /** Checks whether an algebraic number is zero. */
  bool is_zero(const AlgebraicNumber& an);
  /** Checks whether an algebraic number is one. */
  bool is_one(const AlgebraicNumber& an);

  /** Compute the ceiling of an algebraic number. */
  Integer ceil(const AlgebraicNumber& an);
  /** Compute the floor of an algebraic number. */
  Integer floor(const AlgebraicNumber& an);

}  // namespace poly
