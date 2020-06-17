#include "polyxx/dyadic_interval.h"

#include <iostream>

namespace poly {

  DyadicInterval::DyadicInterval() {
    lp_dyadic_interval_construct_zero(get_internal());
  }
  DyadicInterval::DyadicInterval(const DyadicRational& dr) {
    lp_dyadic_interval_construct_point(get_internal(), dr.get_internal());
  }
  DyadicInterval::DyadicInterval(const DyadicRational& a,
                                 const DyadicRational& b)
      : DyadicInterval(a, true, b, true) {}
  DyadicInterval::DyadicInterval(const DyadicRational& a, bool a_open,
                                 const DyadicRational& b, bool b_open) {
    lp_dyadic_interval_construct(get_internal(), a.get_internal(), a_open,
                                 b.get_internal(), b_open);
  }
  DyadicInterval::DyadicInterval(const Integer& i)
      : DyadicInterval(DyadicRational(i)) {}
  DyadicInterval::DyadicInterval(const Integer& a, const Integer& b)
      : DyadicInterval(a, true, b, true) {}
  DyadicInterval::DyadicInterval(const Integer& a, bool a_open,
                                 const Integer& b, bool b_open) {
    lp_dyadic_interval_construct_from_integer(get_internal(), a.get_internal(),
                                              a_open, b.get_internal(), b_open);
  }
  DyadicInterval::DyadicInterval(long a, long b)
      : DyadicInterval(a, true, b, true) {}
  DyadicInterval::DyadicInterval(long a, bool a_open, long b, bool b_open) {
    lp_dyadic_interval_construct_from_int(get_internal(), a, a_open, b, b_open);
  }

  DyadicInterval::DyadicInterval(const DyadicInterval& i) {
    lp_dyadic_interval_construct_copy(get_internal(), i.get_internal());
  }
  DyadicInterval::DyadicInterval(DyadicInterval&& i) {
    lp_dyadic_interval_construct_copy(get_internal(), i.get_internal());
  }

  DyadicInterval::~DyadicInterval() {
    lp_dyadic_interval_destruct(get_internal());
  }

  DyadicInterval& DyadicInterval::operator=(const DyadicInterval& i) {
    lp_dyadic_interval_construct_copy(get_internal(), i.get_internal());
    return *this;
  }
  DyadicInterval& DyadicInterval::operator=(DyadicInterval&& i) {
    lp_dyadic_interval_construct_copy(get_internal(), i.get_internal());
    return *this;
  }

  lp_dyadic_interval_t* DyadicInterval::get_internal() { return &mInterval; }

  const lp_dyadic_interval_t* DyadicInterval::get_internal() const {
    return &mInterval;
  }

  void DyadicInterval::collapse(const DyadicRational& dr) {
    lp_dyadic_interval_collapse_to(get_internal(), dr.get_internal());
  }
  void DyadicInterval::set_lower(const DyadicRational& dr, bool open) {
    lp_dyadic_interval_set_a(get_internal(), dr.get_internal(), open);
  }
  void DyadicInterval::set_upper(const DyadicRational& dr, bool open) {
    lp_dyadic_interval_set_b(get_internal(), dr.get_internal(), open);
  }
  void DyadicInterval::scale(int n) {
    lp_dyadic_interval_scale(get_internal(), n);
  }

  std::ostream& operator<<(std::ostream& os, const DyadicInterval& i) {
    os << (i.get_internal()->a_open ? "( " : "[ ");
    stream_ptr(os, lp_dyadic_rational_to_string(&(i.get_internal()->a)));
    os << " ; ";
    stream_ptr(os, lp_dyadic_rational_to_string(&(i.get_internal()->b)));
    os << (i.get_internal()->b_open ? " )" : " ]");
    return os;
  }

  std::pair<DyadicInterval, DyadicInterval> split(const DyadicInterval& di,
                                                  bool left_open,
                                                  bool right_open) {
    std::pair<DyadicInterval, DyadicInterval> res;
    lp_dyadic_interval_construct_from_split(
        res.first.get_internal(), res.second.get_internal(), di.get_internal(),
        left_open ? 1 : 0, right_open ? 1 : 0);
    return res;
  }

  void swap(DyadicInterval& lhs, DyadicInterval& rhs) {
    lp_dyadic_interval_swap(lhs.get_internal(), rhs.get_internal());
  }

  bool operator==(const DyadicInterval& lhs, const DyadicInterval& rhs) {
    return lp_dyadic_interval_equals(lhs.get_internal(), rhs.get_internal());
  }
  bool operator!=(const DyadicInterval& lhs, const DyadicInterval& rhs) {
    return !(lhs == rhs);
  }

  bool contains(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_contains_dyadic_rational(lhs.get_internal(),
                                                       rhs.get_internal());
  }

  bool contains_zero(const DyadicInterval& lhs) {
    return lp_dyadic_interval_contains_zero(lhs.get_internal());
  }

  bool disjoint(const DyadicInterval& lhs, const DyadicInterval& rhs) {
    return lp_dyadic_interval_disjoint(lhs.get_internal(), rhs.get_internal());
  }

  bool is_point(const DyadicInterval& di) {
    return lp_dyadic_interval_is_point(di.get_internal());
  }
  DyadicRational get_point(const DyadicInterval& di) {
    return DyadicRational(lp_dyadic_interval_get_point(di.get_internal()));
  }
  const DyadicRational& get_lower(const DyadicInterval& i) {
    return *detail::cast_from(&i.get_internal()->a);
  }
  const DyadicRational& get_upper(const DyadicInterval& i) {
    if (is_point(i)) {
      return *detail::cast_from(&i.get_internal()->a);
    } else {
      return *detail::cast_from(&i.get_internal()->b);
    }
  }
  int log_size(const DyadicInterval& di) {
    return lp_dyadic_interval_size(di.get_internal());
  }

  int sgn(const DyadicInterval& di) {
    return lp_dyadic_interval_sgn(di.get_internal());
  }

  bool operator==(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicInterval& lhs, const Integer& rhs) {
    return lp_dyadic_interval_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) >= 0;
  }

  bool operator==(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Integer& lhs, const DyadicInterval& rhs) {
    return rhs <= lhs;
  }

  bool operator==(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicInterval& lhs, const DyadicRational& rhs) {
    return lp_dyadic_interval_cmp_dyadic_rational(lhs.get_internal(),
                                                  rhs.get_internal()) >= 0;
  }

  bool operator==(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs != lhs;
  }
  bool operator<(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const DyadicRational& lhs, const DyadicInterval& rhs) {
    return rhs <= lhs;
  }

  bool operator==(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicInterval& lhs, const Rational& rhs) {
    return lp_dyadic_interval_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) >= 0;
  }

  bool operator==(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Rational& lhs, const DyadicInterval& rhs) {
    return rhs <= lhs;
  }

}  // namespace poly
