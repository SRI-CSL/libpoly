#pragma once

#include <iosfwd>

#include "../interval.h"
#include "value.h"

namespace poly {

  /**
   * Implements a wrapper for lp_interval_t.
   */
  class Interval {
   private:
    /** The actual interval. */
    lp_interval_t mInterval;

   public:
    /** Construct from an internal lp_interval_t pointer. */
    explicit Interval(const lp_interval_t* i);
    /** Construct zero point interval. */
    Interval();
    /** Construct an open interval from the given two values and bound types. */
    Interval(const Value& a, bool a_open, const Value& b, bool b_open);
    /** Construct an open interval from the given two values. */
    Interval(const Value& a, const Value& b);
    /** Construct a point interval from the given value. */
    Interval(const Value& a);
    /** Copy from the given Interval. */
    Interval(const Interval& i);
    /** Move from the given Interval. */
    Interval(Interval&& i);
    /** Custom destructor. */
    ~Interval();
    /** Copy from the given Interval. */
    Interval& operator=(const Interval& i);
    /** Move from the given Interval. */
    Interval& operator=(Interval&& i);

    /** Get a non-const pointer to the internal lp_interval_t. Handle with
     * care! */
    lp_interval_t* get_internal();
    /** Get a const pointer to the internal lp_interval_t. */
    const lp_interval_t* get_internal() const;

    /** Collapse this interval to a single point. */
    void collapse(const Value& v);
    /** The the lower bound. */
    void set_lower(const Value& v, bool open);
    /** The the upper bound. */
    void set_upper(const Value& v, bool open);

    /** Construct the full interval. */
    static Interval full();
  };

  /** Make sure that we can cast between Interval and lp_interval_t. */
  static_assert(sizeof(Interval) == sizeof(lp_interval_t),
                "Please check the size of Interval.");
  namespace detail {
    /** Non-const cast from an Interval to a lp_interval_t. */
    inline lp_interval_t* cast_to(Interval* i) {
      return reinterpret_cast<lp_interval_t*>(i);
    }
    /** Const cast from an Interval to a lp_interval_t. */
    inline const lp_interval_t* cast_to(const Interval* i) {
      return reinterpret_cast<const lp_interval_t*>(i);
    }
    /** Non-const cast from a lp_interval_t to an Interval. */
    inline Interval* cast_from(lp_interval_t* i) {
      return reinterpret_cast<Interval*>(i);
    }
    /** Const cast from a lp_interval_t to an Interval. */
    inline const Interval* cast_from(const lp_interval_t* i) {
      return reinterpret_cast<const Interval*>(i);
    }
  }  // namespace detail

  /** Swap two intervals. */
  void swap(Interval& lhs, Interval& rhs);

  /** Check whether an interval contains some value. */
  bool contains(const Interval& i, const Value& v);
  /** Get an approximation of the log interval size. */
  int log_size(const Interval& i);

  /** Stream the given Interval to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Interval& i);

  /** Check whether an interval is a point. */
  bool is_point(const Interval& i);
  /** Check whether an interval is the full interval. */
  bool is_full(const Interval& i);
  /** Get the value of a point interval. Assumes is_point(i). */
  const Value& get_point(const Interval& i);
  /** Get the lower bound of an interval. */
  const Value& get_lower(const Interval& i);
  /** Get the upper bound of an interval. */
  const Value& get_upper(const Interval& i);

  /** Pick some value from an interval. */
  Value pick_value(const Interval& i);

  /** Compare the lower bounds of two intervals. */
  int compare_lower(const Interval& lhs, const Interval& rhs);
  /** Compare the upper bounds of two intervals. */
  int compare_upper(const Interval& lhs, const Interval& rhs);

  /** Compare two intervals. */
  bool operator==(const Interval& lhs, const Interval& rhs);
  /** Compare two intervals. */
  bool operator!=(const Interval& lhs, const Interval& rhs);
  /** Compare two intervals. */
  bool operator<(const Interval& lhs, const Interval& rhs);
  /** Compare two intervals. */
  bool operator<=(const Interval& lhs, const Interval& rhs);
  /** Compare two intervals. */
  bool operator>(const Interval& lhs, const Interval& rhs);
  /** Compare two intervals. */
  bool operator>=(const Interval& lhs, const Interval& rhs);

  /** Compute the sign of an interval. */
  int sgn(const Interval& i);
  /** Compute i^n. */
  Interval pow(const Interval& i, unsigned n);
  /** Add two intervals. */
  Interval operator+(const Interval& lhs, const Interval& rhs);
  /** Multiply two intervals. */
  Interval operator*(const Interval& lhs, const Interval& rhs);

}  // namespace poly
