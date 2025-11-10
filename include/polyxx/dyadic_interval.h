#pragma once

#include <iosfwd>

#include "../dyadic_interval.h"
#include "dyadic_rational.h"
#include "integer.h"

namespace poly {

  /**
   * Implements a wrapper for lp_dyadic_interval_t.
   */
  class DyadicInterval {
   private:
    /** The actual interval. */
    lp_dyadic_interval_t mInterval;

   public:
    /** Construct from an internal lp_dyadic_interval_t pointer. */
    explicit DyadicInterval(const lp_dyadic_interval_t* di);
    /** Construct a zero interval [0;0]. */
    DyadicInterval();
    /** Construct a point interval. */
    DyadicInterval(const DyadicRational& dr);
    /** Construct an open interval. */
    DyadicInterval(const DyadicRational& a, const DyadicRational& b);
    /** Construct an interval from the given bounds. */
    DyadicInterval(const DyadicRational& a, bool a_open,
                   const DyadicRational& b, bool b_open);
    /** Construct a point interval. */
    DyadicInterval(const Integer& i);
    /** Construct an open interval. */
    DyadicInterval(const Integer& a, const Integer& b);
    /** Construct an interval from the given bounds. */
    DyadicInterval(const Integer& a, bool a_open,
                   const Integer& b, bool b_open);
    /** Construct a point interval. */
    explicit DyadicInterval(long i);
    /** Construct an open interval. */
    explicit DyadicInterval(long a, long b);
    /** Construct an interval from the given bounds. */
    explicit DyadicInterval(long a, bool a_open, long b, bool b_open);
    /** Copy from a DyadicInterval. */
    DyadicInterval(const DyadicInterval& i);
    /** Move from a DyadicInterval. */
    DyadicInterval(DyadicInterval&& i);
    /** Custom destructor. */
    ~DyadicInterval();
    /** Copy from a DyadicInterval. */
    DyadicInterval& operator=(const DyadicInterval& i);
    /** Move from a DyadicInterval. */
    DyadicInterval& operator=(DyadicInterval&& i);

    /** Get a non-const pointer to the internal lp_dyadic_interval_t. Handle
     * with care! */
    lp_dyadic_interval_t* get_internal();
    /** Get a const pointer to the internal lp_dyadic_interval_t. */
    const lp_dyadic_interval_t* get_internal() const;

    /** Collapse this interval to a single point. */
    void collapse(const DyadicRational& dr);
    /** The lower bound. */
    void set_lower(const DyadicRational& dr, bool open);
    /** The upper bound. */
    void set_upper(const DyadicRational& dr, bool open);
    /** The interval by 2^n. */
    void scale(int n);
  };

  /** Make sure that we can cast between DyadicInterval and
   * lp_dyadic_interval_t. */
  static_assert(sizeof(DyadicInterval) == sizeof(lp_dyadic_interval_t),
                "Please check the size of DyadicInterval.");
  namespace detail {
    /** Non-const cast from an DyadicInterval to a lp_dyadic_interval_t. */
    inline lp_dyadic_interval_t* cast_to(DyadicInterval* i) {
      return reinterpret_cast<lp_dyadic_interval_t*>(i);
    }
    /** Const cast from an DyadicInterval to a lp_dyadic_interval_t. */
    inline const lp_dyadic_interval_t* cast_to(const DyadicInterval* i) {
      return reinterpret_cast<const lp_dyadic_interval_t*>(i);
    }
    /** Non-const cast from a lp_dyadic_interval_t to an DyadicInterval. */
    inline DyadicInterval* cast_from(lp_dyadic_interval_t* i) {
      return reinterpret_cast<DyadicInterval*>(i);
    }
    /** Const cast from a lp_dyadic_interval_t to an DyadicInterval. */
    inline const DyadicInterval* cast_from(const lp_dyadic_interval_t* i) {
      return reinterpret_cast<const DyadicInterval*>(i);
    }
  }  // namespace detail

  /** Stream the given DyadicInterval to an output stream. */
  std::ostream& operator<<(std::ostream& os, const DyadicInterval& i);

  /** Split the given interval. */
  std::pair<DyadicInterval, DyadicInterval> split(const DyadicInterval& di,
                                                  bool left_open,
                                                  bool right_open);

  /** Swap two intervals. */
  void swap(DyadicInterval& lhs, DyadicInterval& rhs);

  /** Compare two intervals. */
  bool operator==(const DyadicInterval& lhs, const DyadicInterval& rhs);
  /** Compare two intervals. */
  bool operator!=(const DyadicInterval& lhs, const DyadicInterval& rhs);

  /** Check whether the interval contains a given point. */
  bool contains(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Check whether the interval contains zero. */
  bool contains_zero(const DyadicInterval& lhs);

  /** Check whether two intervals are disjoint. */
  bool disjoint(const DyadicInterval& lhs, const DyadicInterval& rhs);

  /** Check whether an interval is a point. */
  bool is_point(const DyadicInterval& di);
  /** Get the point from a point interval. */
  DyadicRational get_point(const DyadicInterval& di);
  /** Get the lower bound. */
  const DyadicRational& get_lower(const DyadicInterval& i);
  /** Get the upper bound. */
  const DyadicRational& get_upper(const DyadicInterval& i);

  /** Get approximate log size. */
  int log_size(const DyadicInterval& di);

  /** Get sign of the interval. */
  int sgn(const DyadicInterval& di);

  /** Compare an interval and an integer. */
  bool operator==(const DyadicInterval& lhs, const Integer& rhs);
  /** Compare an interval and an integer. */
  bool operator!=(const DyadicInterval& lhs, const Integer& rhs);
  /** Compare an interval and an integer. */
  bool operator<(const DyadicInterval& lhs, const Integer& rhs);
  /** Compare an interval and an integer. */
  bool operator<=(const DyadicInterval& lhs, const Integer& rhs);
  /** Compare an interval and an integer. */
  bool operator>(const DyadicInterval& lhs, const Integer& rhs);
  /** Compare an interval and an integer. */
  bool operator>=(const DyadicInterval& lhs, const Integer& rhs);

  /** Compare an integer and an interval. */
  bool operator==(const Integer& lhs, const DyadicInterval& rhs);
  /** Compare an integer and an interval. */
  bool operator!=(const Integer& lhs, const DyadicInterval& rhs);
  /** Compare an integer and an interval. */
  bool operator<(const Integer& lhs, const DyadicInterval& rhs);
  /** Compare an integer and an interval. */
  bool operator<=(const Integer& lhs, const DyadicInterval& rhs);
  /** Compare an integer and an interval. */
  bool operator>(const Integer& lhs, const DyadicInterval& rhs);
  /** Compare an integer and an interval. */
  bool operator>=(const Integer& lhs, const DyadicInterval& rhs);

  /** Compare an interval and a rational. */
  bool operator==(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Compare an interval and a rational. */
  bool operator!=(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Compare an interval and a rational. */
  bool operator<(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Compare an interval and a rational. */
  bool operator<=(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Compare an interval and a rational. */
  bool operator>(const DyadicInterval& lhs, const DyadicRational& rhs);
  /** Compare an interval and a rational. */
  bool operator>=(const DyadicInterval& lhs, const DyadicRational& rhs);

  /** Compare a rational and an interval. */
  bool operator==(const DyadicRational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator!=(const DyadicRational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator<(const DyadicRational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator<=(const DyadicRational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator>(const DyadicRational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator>=(const DyadicRational& lhs, const DyadicInterval& rhs);

  /** Compare an interval and a rational. */
  bool operator==(const DyadicInterval& lhs, const Rational& rhs);
  /** Compare an interval and a rational. */
  bool operator!=(const DyadicInterval& lhs, const Rational& rhs);
  /** Compare an interval and a rational. */
  bool operator<(const DyadicInterval& lhs, const Rational& rhs);
  /** Compare an interval and a rational. */
  bool operator<=(const DyadicInterval& lhs, const Rational& rhs);
  /** Compare an interval and a rational. */
  bool operator>(const DyadicInterval& lhs, const Rational& rhs);
  /** Compare an interval and a rational. */
  bool operator>=(const DyadicInterval& lhs, const Rational& rhs);

  /** Compare a rational and an interval. */
  bool operator==(const Rational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator!=(const Rational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator<(const Rational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator<=(const Rational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator>(const Rational& lhs, const DyadicInterval& rhs);
  /** Compare a rational and an interval. */
  bool operator>=(const Rational& lhs, const DyadicInterval& rhs);

}  // namespace poly
