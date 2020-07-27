#include "polyxx/interval.h"

#include <cassert>

namespace poly {

  Interval::Interval(const lp_interval_t* i) {
    lp_interval_construct_copy(get_internal(), i);
  }

  Interval::Interval() { lp_interval_construct_zero(get_internal()); }

  Interval::Interval(const Value& a, bool a_open, const Value& b, bool b_open) {
    lp_interval_construct(get_internal(), a.get_internal(), a_open ? 1 : 0,
                          b.get_internal(), b_open ? 1 : 0);
  }

  Interval::Interval(const Value& a, const Value& b)
      : Interval(a, true, b, true) {}

  Interval::Interval(const Value& a) {
    lp_interval_construct_point(get_internal(), a.get_internal());
  }

  Interval::Interval(const Interval& i) {
    lp_interval_construct_copy(get_internal(), i.get_internal());
  }
  Interval::Interval(Interval&& i) {
    lp_interval_construct_copy(get_internal(), i.get_internal());
  }

  Interval::~Interval() { lp_interval_destruct(get_internal()); }

  Interval& Interval::operator=(const Interval& i) {
    lp_interval_construct_copy(get_internal(), i.get_internal());
    return *this;
  }
  Interval& Interval::operator=(Interval&& i) {
    swap(*this, i);
    return *this;
  }

  lp_interval_t* Interval::get_internal() { return &mInterval; }

  const lp_interval_t* Interval::get_internal() const { return &mInterval; }

  void Interval::collapse(const Value& v) {
    lp_interval_collapse_to(get_internal(), v.get_internal());
  }
  void Interval::set_lower(const Value& v, bool open) {
    lp_interval_set_a(get_internal(), v.get_internal(), open ? 1 : 0);
  }
  void Interval::set_upper(const Value& v, bool open) {
    lp_interval_set_b(get_internal(), v.get_internal(), open ? 1 : 0);
  }

  Interval Interval::full() {
    Interval res;
    lp_interval_construct_full(res.get_internal());
    return res;
  }

  void swap(Interval& lhs, Interval& rhs) {
    lp_interval_swap(lhs.get_internal(), rhs.get_internal());
  }

  bool contains(const Interval& i, const Value& v) {
    return lp_interval_contains(i.get_internal(), v.get_internal());
  }
  int log_size(const Interval& i) {
    return lp_interval_size_approx(i.get_internal());
  }

  std::ostream& operator<<(std::ostream& os, const Interval& i) {
    return stream_ptr(os, lp_interval_to_string(i.get_internal()));
  }

  bool is_point(const Interval& i) {
    return lp_interval_is_point(i.get_internal());
  }
  bool is_full(const Interval& i) {
    return lp_interval_is_full(i.get_internal());
  }
  const Value& get_point(const Interval& i) {
    return *detail::cast_from(lp_interval_get_point(i.get_internal()));
  }
  const Value& get_lower(const Interval& i) {
    return *detail::cast_from(lp_interval_get_lower_bound(i.get_internal()));
  }
  bool get_lower_open(const Interval& i) {
    return !i.get_internal()->is_point && i.get_internal()->a_open;
  }
  const Value& get_upper(const Interval& i) {
    return *detail::cast_from(lp_interval_get_upper_bound(i.get_internal()));
  }
  bool get_upper_open(const Interval& i) {
    return !i.get_internal()->is_point && i.get_internal()->b_open;
  }

  Value pick_value(const Interval& i) {
    Value res;
    lp_interval_pick_value(i.get_internal(), res.get_internal());
    return res;
  }

  int compare_lower(const Interval& lhs, const Interval& rhs) {
    return lp_interval_cmp_lower_bounds(lhs.get_internal(), rhs.get_internal());
  }
  int compare_upper(const Interval& lhs, const Interval& rhs) {
    return lp_interval_cmp_upper_bounds(lhs.get_internal(), rhs.get_internal());
  }

  lp_interval_cmp_t compare(const Interval& lhs, const Interval& rhs) {
    return lp_interval_cmp(lhs.get_internal(), rhs.get_internal());
  }

  bool operator==(const Interval& lhs, const Interval& rhs) {
    lp_interval_cmp_t res = compare(lhs, rhs);
    return res == LP_INTERVAL_CMP_EQ;
  }
  bool operator!=(const Interval& lhs, const Interval& rhs) {
    return !(lhs == rhs);
  }
  bool operator<(const Interval& lhs, const Interval& rhs) {
    int l = compare_lower(lhs, rhs);
    if (l != 0) return l < 0;
    return compare_upper(lhs, rhs) < 0;
  }
  bool operator<=(const Interval& lhs, const Interval& rhs) {
    int l = compare_lower(lhs, rhs);
    if (l != 0) return l < 0;
    return compare_upper(lhs, rhs) <= 0;
  }
  bool operator>(const Interval& lhs, const Interval& rhs) { return rhs < lhs; }
  bool operator>=(const Interval& lhs, const Interval& rhs) {
    return rhs <= lhs;
  }

  int sgn(const Interval& i) { return lp_interval_sgn(i.get_internal()); }
  Interval pow(const Interval& i, unsigned n) {
    Interval res;
    lp_interval_pow(res.get_internal(), i.get_internal(), n);
    return res;
  }
  Interval operator+(const Interval& lhs, const Interval& rhs) {
    Interval res;
    lp_interval_add(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }
  Interval operator*(const Interval& lhs, const Interval& rhs) {
    Interval res;
    lp_interval_mul(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }

}  // namespace poly
