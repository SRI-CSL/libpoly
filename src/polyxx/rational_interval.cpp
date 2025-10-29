#include "polyxx/rational_interval.h"

#include <cassert>
#include <iostream>

namespace poly {

  RationalInterval::RationalInterval(const lp_rational_interval_t* i) {
    lp_rational_interval_construct_copy(get_internal(), i);
  }

  RationalInterval::RationalInterval() : RationalInterval(Rational(0)) {}

  RationalInterval::RationalInterval(const Rational& a, bool a_open,
                                     const Rational& b, bool b_open) {
    lp_rational_interval_construct(get_internal(), a.get_internal(),
                                   a_open ? 1 : 0, b.get_internal(),
                                   b_open ? 1 : 0);
  }

  RationalInterval::RationalInterval(const Rational& a, const Rational& b)
      : RationalInterval(a, true, b, true) {}

  RationalInterval::RationalInterval(const Rational& a) {
    lp_rational_interval_construct_point(get_internal(), a.get_internal());
  }

  RationalInterval::RationalInterval(const RationalInterval& i) {
    lp_rational_interval_construct_copy(get_internal(), i.get_internal());
  }

  RationalInterval::RationalInterval(const DyadicRational& a, bool a_open,
                                     const DyadicRational& b, bool b_open) {
    lp_rational_interval_construct_from_dyadic(get_internal(), a.get_internal(),
                                               a_open ? 1 : 0, b.get_internal(),
                                               b_open ? 1 : 0);
  }
  RationalInterval::RationalInterval(const DyadicRational& a,
                                     const DyadicRational& b)
      : RationalInterval(a, true, b, true) {}
  RationalInterval::RationalInterval(const DyadicInterval& i) {
    lp_rational_interval_construct_from_dyadic_interval(get_internal(),
                                                        i.get_internal());
  }

  RationalInterval::~RationalInterval() {
    lp_rational_interval_destruct(&mInterval);
  }

  RationalInterval& RationalInterval::operator=(const RationalInterval& i) {
    lp_rational_interval_construct_copy(get_internal(), i.get_internal());
    return *this;
  }
  RationalInterval& RationalInterval::operator=(RationalInterval&& i) {
    swap(*this, i);
    return *this;
  }

  lp_rational_interval_t* RationalInterval::get_internal() {
    return &mInterval;
  }

  const lp_rational_interval_t* RationalInterval::get_internal() const {
    return &mInterval;
  }

  void swap(RationalInterval& lhs, RationalInterval& rhs) {
    lp_rational_interval_swap(lhs.get_internal(), rhs.get_internal());
  }

  std::ostream& operator<<(std::ostream& os, const RationalInterval& i) {
    if (i.get_internal()->is_point) {
      assert(!i.get_internal()->a_open && !i.get_internal()->b_open);
      os << "[ ";
      stream_ptr(os, lp_rational_to_string(&(i.get_internal()->a)));
      os << " ; ";
      stream_ptr(os, lp_rational_to_string(&(i.get_internal()->a)));
      return os << " ]";
    }
    os << (i.get_internal()->a_open ? "( " : "[ ");
    stream_ptr(os, lp_rational_to_string(&(i.get_internal()->a)));
    os << " ; ";
    stream_ptr(os, lp_rational_to_string(&(i.get_internal()->b)));
    os << (i.get_internal()->b_open ? " )" : " ]");
    return os;
  }

  bool contains(const RationalInterval& ri, const AlgebraicNumber& an) {
    assert(false && "Not yet implemented in the library.");
    return lp_rational_interval_contains_algebraic_number(ri.get_internal(),
                                                          an.get_internal());
  }
  bool contains(const RationalInterval& ri, const DyadicRational& dr) {
    assert(false && "Not yet implemented in the library.");
    return lp_rational_interval_contains_dyadic_rational(ri.get_internal(),
                                                         dr.get_internal());
  }
  bool contains(const RationalInterval& ri, const Integer& i) {
    assert(false && "Not yet implemented in the library.");
    return lp_rational_interval_contains_integer(ri.get_internal(),
                                                 i.get_internal());
  }
  bool contains(const RationalInterval& ri, const Rational& r) {
    return lp_rational_interval_contains_rational(ri.get_internal(),
                                                  r.get_internal());
  }
  bool contains(const RationalInterval& ri, const Value& v) {
    return lp_rational_interval_contains_value(ri.get_internal(),
                                               v.get_internal());
  }
  bool contains_zero(const RationalInterval& ri) {
    return lp_rational_interval_contains_zero(ri.get_internal());
  }

  bool is_point(const RationalInterval& ri) {
    return lp_rational_interval_is_point(ri.get_internal());
  }
  const Rational& get_point(const RationalInterval& i) {
    return *detail::cast_from(lp_rational_interval_get_point(i.get_internal()));
  }
  const Rational& get_lower(const RationalInterval& ri) {
    return *detail::cast_from(&ri.get_internal()->a);
  }
  const Rational& get_upper(const RationalInterval& ri) {
    if (is_point(ri)) {
      return get_lower(ri);
    }
    return *detail::cast_from(&ri.get_internal()->b);
  }

  int sgn(const RationalInterval& ri) {
    return lp_rational_interval_sgn(ri.get_internal());
  }

}  // namespace poly
