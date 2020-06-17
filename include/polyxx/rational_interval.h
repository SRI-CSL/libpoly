#pragma once

#include <iosfwd>

#include "../rational_interval.h"
#include "algebraic_number.h"
#include "dyadic_interval.h"
#include "integer.h"
#include "rational.h"
#include "value.h"

namespace poly {

  /**
   * Implements a wrapper for lp_rational_interval_t.
   */
  class RationalInterval {
   private:
    /** The actual interval. */
    lp_rational_interval_t mInterval;

   public:
    /** Construct from an internal lp_rational_interval_t pointer. */
    explicit RationalInterval(const lp_rational_interval_t* ri);
    /** Construct zero point interval. */
    RationalInterval();
    /** Construct from the given bounds. */
    RationalInterval(const Rational& a, bool a_open, const Rational& b,
                     bool b_open);
    /** Construct from the given bounds. */
    RationalInterval(const Rational& a, const Rational& b);
    /** Construct point inverval. */
    RationalInterval(const Rational& a);
    /** Construct from the given bounds. */
    RationalInterval(const DyadicRational& a, bool a_open,
                     const DyadicRational& b, bool b_open);
    /** Construct from the given bounds. */
    RationalInterval(const DyadicRational& a, const DyadicRational& b);
    /** Copy from the DyadicInterval. */
    RationalInterval(const DyadicInterval& i);
    /** Copy from the RationalInterval. */
    RationalInterval(const RationalInterval& i);
    /** Move from the RationalInterval. */
    RationalInterval(RationalInterval&& i);
    /** Custom destructor. */
    ~RationalInterval();

    /** Copy from the RationalInterval. */
    RationalInterval& operator=(const RationalInterval& i);
    /** Move from the RationalInterval. */
    RationalInterval& operator=(RationalInterval&& i);

    /** Get a non-const pointer to the internal lp_interval_t. Handle with
     * care! */
    lp_rational_interval_t* get_internal();
    /** Get a const pointer to the internal lp_interval_t. */
    const lp_rational_interval_t* get_internal() const;
  };

  /** Swap two intervals. */
  void swap(RationalInterval& lhs, RationalInterval& rhs);

  /** Stream the given Interval to an output stream. */
  std::ostream& operator<<(std::ostream& os, const RationalInterval& i);

  /** Checks whether an interval contains an algebraic number. */
  bool contains(const RationalInterval& ri, const AlgebraicNumber& an);
  /** Checks whether an interval contains a dyadic rational. */
  bool contains(const RationalInterval& ri, const DyadicRational& dr);
  /** Checks whether an interval contains an integer. */
  bool contains(const RationalInterval& ri, const Integer& i);
  /** Checks whether an interval contains a rational. */
  bool contains(const RationalInterval& ri, const Rational& r);
  /** Checks whether an interval contains a value. */
  bool contains(const RationalInterval& ri, const Value& v);
  /** Checks whether an interval contains zero. */
  bool contains_zero(const RationalInterval& ri);

  /** Checks whether an interval is a point interval. */
  bool is_point(const RationalInterval& ri);
  /** Returns the value of a point interval. */
  const Rational& get_point(const RationalInterval& ri);
  /** Returns the lower bound of an interval. */
  const Rational& get_lower(const RationalInterval& ri);
  /** Returns the upper bound of an interval. */
  const Rational& get_upper(const RationalInterval& ri);

  /** Compute the sign of an interval. */
  int sgn(const RationalInterval& ri);

}  // namespace poly
