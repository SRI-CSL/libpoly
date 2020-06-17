#pragma once

#include <iosfwd>

#include "../value.h"
#include "algebraic_number.h"
#include "dyadic_rational.h"
#include "integer.h"
#include "rational.h"

namespace poly {

  /**
   * Implements a wrapper for lp_value_t.
   * A Value can hold an Integer, a DyadicRational, a Rational, an
   * AlgebraicNumber or a special value for nothing, plus infinity or minus
   * infinity.
   */
  class Value {
    /** The actual value. */
    lp_value_t mValue;

    /** Internal wrapper of the general Value constructor. */
    Value(lp_value_type_t type, const void* data);

   public:
    /** Create from a lp_value_t, creating a copy. */
    explicit Value(const lp_value_t* val);
    /** Construct a none value. */
    Value();
    /** Construct from a native integer. */
    Value(long i);
    /** Copy from the given Value. */
    Value(const Value& val);
    /** Move from the given Value. */
    Value(Value&& val);

    /** Construct from an algebraic number. */
    Value(const AlgebraicNumber& an);
    /** Construct from a dyadic rational. */
    Value(const DyadicRational& dr);
    /** Construct from an integer. */
    Value(const Integer& i);
    /** Construct from a rational. */
    Value(const Rational& r);

    /** Custom destructor. */
    ~Value();

    /** Copy from the given Value. */
    Value& operator=(const Value& v);
    /** Move from the given Value. */
    Value& operator=(Value&& v);

    /** Get a non-const pointer to the internal lp_value_t. Handle with care! */
    lp_value_t* get_internal();
    /** Get a const pointer to the internal lp_value_t. */
    const lp_value_t* get_internal() const;

    /** Return -infty */
    static Value minus_infty();
    /** Return +infty */
    static Value plus_infty();
  };

  /** Make sure that we can cast between Value and lp_value_t. */
  static_assert(sizeof(Value) == sizeof(lp_value_t),
                "Please check the size of Value.");
  namespace detail {
    /** Non-const cast from a Value to a lp_value_t. */
    inline lp_value_t* cast_to(Value* i) {
      return reinterpret_cast<lp_value_t*>(i);
    }
    /** Const cast from a Value to a lp_value_t. */
    inline const lp_value_t* cast_to(const Value* i) {
      return reinterpret_cast<const lp_value_t*>(i);
    }
    /** Non-const cast from a lp_value_t to a Value. */
    inline Value* cast_from(lp_value_t* i) {
      return reinterpret_cast<Value*>(i);
    }
    /** Const cast from a lp_value_t to a Value. */
    inline const Value* cast_from(const lp_value_t* i) {
      return reinterpret_cast<const Value*>(i);
    }
  }  // namespace detail

  /** Compare two values. */
  bool operator==(const Value& lhs, const Value& rhs);
  /** Compare two values. */
  bool operator!=(const Value& lhs, const Value& rhs);
  /** Compare two values. */
  bool operator<(const Value& lhs, const Value& rhs);
  /** Compare two values. */
  bool operator<=(const Value& lhs, const Value& rhs);
  /** Compare two values. */
  bool operator>(const Value& lhs, const Value& rhs);
  /** Compare two values. */
  bool operator>=(const Value& lhs, const Value& rhs);

  /** Compare a value and a rational. */
  bool operator==(const Value& lhs, const Rational& rhs);
  /** Compare a value and a rational. */
  bool operator!=(const Value& lhs, const Rational& rhs);
  /** Compare a value and a rational. */
  bool operator<(const Value& lhs, const Rational& rhs);
  /** Compare a value and a rational. */
  bool operator<=(const Value& lhs, const Rational& rhs);
  /** Compare a value and a rational. */
  bool operator>(const Value& lhs, const Rational& rhs);
  /** Compare a value and a rational. */
  bool operator>=(const Value& lhs, const Rational& rhs);

  /** Compare a rational and a value. */
  bool operator==(const Rational& lhs, const Value& rhs);
  /** Compare a rational and a value. */
  bool operator!=(const Rational& lhs, const Value& rhs);
  /** Compare a rational and a value. */
  bool operator<(const Rational& lhs, const Value& rhs);
  /** Compare a rational and a value. */
  bool operator<=(const Rational& lhs, const Value& rhs);
  /** Compare a rational and a value. */
  bool operator>(const Rational& lhs, const Value& rhs);
  /** Compare a rational and a value. */
  bool operator>=(const Rational& lhs, const Value& rhs);

  /** Stream the given Value to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Value& v);

  /** Swap two values. */
  void swap(Value& lhs, Value& rhs);
  /** Compute the hash of a value.
   * To obtain a consistent hash among all representations, the floor of the
   * stored number is hashed.
   */
  std::size_t hash(const Value& v);

  /** Returns the sign of a value. */
  int sgn(const Value& v);

  /** Checks whether a value holds an algebraic number. */
  bool is_algebraic_number(const Value& v);
  /** Checks whether a value holds a dyadic rational. */
  bool is_dyadic_rational(const Value& v);
  /** Checks whether a value holds an integer. */
  bool is_integer(const Value& v);
  /** Checks whether a value holds minus infinity. */
  bool is_minus_infinity(const Value& v);
  /** Checks whether a value holds nothing. */
  bool is_none(const Value& v);
  /** Checks whether a value holds plus infinity. */
  bool is_plus_infinity(const Value& v);
  /** Checks whether a value holds a rational. */
  bool is_rational(const Value& v);

  /** Return the algebraic number a value holds. Assumes is_algebraic_number(v).
   */
  const AlgebraicNumber& as_algebraic_number(const Value& v);
  /** Return the dyadic rational a value holds. Assumes is_dyadic_rational(v).
   */
  const DyadicRational& as_dyadic_rational(const Value& v);
  /** Return the integer a value holds. Assumes is_integer(v). */
  const Integer& as_integer(const Value& v);
  /** Return the rational a value holds. Assumes is_rational(v). */
  const Rational& as_rational(const Value& v);

  /** Checks whether the value can be interpreted as an integer.
   * Returns true if v holds an Integer, a Rational or DyadicRational with
   * denominator 1, or an AlgebraicNumber that happens to be integral.
   */
  bool represents_integer(const Value& v);
  /** Checks whether the value can be interpreted as an rational.
   * Returns true if v holds an Integer, a Rational or a DyadicRational, or an
   * AlgebraicNumber that happens to be rational.
   */
  bool represents_rational(const Value& v);
  /** Return the integer represented by a value. Assumes represents_integer(v).
   */
  Integer get_integer(const Value& v);
  /** Return the rational represented by a value. Assumes
   * represents_rational(v). */
  Rational get_rational(const Value& v);
  /** Return a double approximation of a value. */
  double to_double(const Value& v);

  /** Computes the numerator of a value that holds an integer, rational or
   * dyadic rational. */
  Integer numerator(const Value& v);
  /** Computes the denominator of a value that holds an integer, rational or
   * dyadic rational. */
  Integer denominator(const Value& v);

  /** Computes the ceiling of a value. */
  Integer ceil(const Value& v);
  /** Computes the floor of a value. */
  Integer floor(const Value& v);

  /** Computs some value within the given interval. */
  Value value_between(const lp_value_t* lhs, bool l_strict,
                      const lp_value_t* rhs, bool r_strict);

  /** Computs some value within the given interval. */
  Value value_between(const Value& lhs, bool l_strict, const Value& rhs,
                      bool r_strict);

  /** Returns the approximate size of the distance of two values. */
  int approximate_size(const Value& lower, const Value& upper);

}  // namespace poly
