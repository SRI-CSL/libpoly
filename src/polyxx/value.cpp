#include "polyxx/value.h"

#include <cassert>

namespace poly {

  Value::Value(lp_value_type_t type, const void* data) {
    lp_value_construct(get_internal(), type, data);
  }

  Value::Value(const lp_value_t* val) {
    lp_value_construct_copy(get_internal(), val);
  }

  Value::Value() { lp_value_construct_none(&mValue); }
  Value::Value(long i) { lp_value_construct_int(get_internal(), i); }
  Value::Value(const Value& val) {
    lp_value_construct_copy(get_internal(), val.get_internal());
  }
  Value::Value(Value&& val) {
    lp_value_construct_copy(get_internal(), val.get_internal());
  }

  Value::Value(const AlgebraicNumber& an)
      : Value(lp_value_type_t::LP_VALUE_ALGEBRAIC, an.get_internal()) {}
  Value::Value(const DyadicRational& dr)
      : Value(lp_value_type_t::LP_VALUE_DYADIC_RATIONAL, dr.get_internal()) {}
  Value::Value(const Integer& i)
      : Value(lp_value_type_t::LP_VALUE_INTEGER, i.get_internal()) {}
  Value::Value(const Rational& r)
      : Value(lp_value_type_t::LP_VALUE_RATIONAL, r.get_internal()) {}

  Value::~Value() {
    lp_value_destruct(get_internal());
  }

  Value& Value::operator=(const Value& v) {
    lp_value_assign(get_internal(), v.get_internal());
    return *this;
  }
  Value& Value::operator=(Value&& v) {
    lp_value_assign(get_internal(), v.get_internal());
    return *this;
  }

  lp_value_t* Value::get_internal() { return &mValue; }
  const lp_value_t* Value::get_internal() const { return &mValue; }

  Value Value::minus_infty() { return Value(LP_VALUE_MINUS_INFINITY, nullptr); }
  Value Value::plus_infty() { return Value(LP_VALUE_PLUS_INFINITY, nullptr); }

  bool operator==(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const Value& lhs, const Value& rhs) {
    return lp_value_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const Value& lhs, const Rational& rhs) {
    return lp_value_cmp_rational(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) == 0;
  }
  bool operator!=(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) != 0;
  }
  bool operator<(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) > 0;
  }
  bool operator<=(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) >= 0;
  }
  bool operator>(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) < 0;
  }
  bool operator>=(const Rational& lhs, const Value& rhs) {
    return lp_value_cmp_rational(rhs.get_internal(), lhs.get_internal()) <= 0;
  }

  std::ostream& operator<<(std::ostream& os, const Value& v) {
    return stream_ptr(os, lp_value_to_string(v.get_internal()));
  }

  void swap(Value& lhs, Value& rhs) {
    lp_value_swap(lhs.get_internal(), rhs.get_internal());
  }
  std::size_t hash(const Value& v) { return lp_value_hash(v.get_internal()); }

  int sgn(const Value& v) { return lp_value_sgn(v.get_internal()); }

  bool is_algebraic_number(const Value& v) {
    return v.get_internal()->type == LP_VALUE_ALGEBRAIC;
  }
  bool is_dyadic_rational(const Value& v) {
    return v.get_internal()->type == LP_VALUE_DYADIC_RATIONAL;
  }
  bool is_integer(const Value& v) {
    return v.get_internal()->type == LP_VALUE_INTEGER;
  }
  bool is_minus_infinity(const Value& v) {
    return v.get_internal()->type == LP_VALUE_MINUS_INFINITY;
  }
  bool is_none(const Value& v) {
    return v.get_internal()->type == LP_VALUE_NONE;
  }
  bool is_plus_infinity(const Value& v) {
    return v.get_internal()->type == LP_VALUE_PLUS_INFINITY;
  }
  bool is_rational(const Value& v) {
    return v.get_internal()->type == LP_VALUE_RATIONAL;
  }

  const AlgebraicNumber& as_algebraic_number(const Value& v) {
    return *detail::cast_from(&v.get_internal()->value.a);
  }
  const DyadicRational& as_dyadic_rational(const Value& v) {
    return *detail::cast_from(&v.get_internal()->value.dy_q);
  }
  const Integer& as_integer(const Value& v) {
    return *detail::cast_from(&v.get_internal()->value.z);
  }
  const Rational& as_rational(const Value& v) {
    return *detail::cast_from(&v.get_internal()->value.q);
  }

  bool represents_integer(const Value& v) {
    return lp_value_is_integer(v.get_internal());
  }
  bool represents_rational(const Value& v) {
    return lp_value_is_rational(v.get_internal());
  }
  Integer get_integer(const Value& v) {
    assert(represents_integer(v));
    Rational res;
    lp_value_get_rational(v.get_internal(), res.get_internal());
    assert(denominator(res) == Rational(1));
    return numerator(res);
  }
  Rational get_rational(const Value& v) {
    assert(represents_rational(v));
    Rational res;
    lp_value_get_rational(v.get_internal(), res.get_internal());
    return res;
  }
  double to_double(const Value& v) {
    return lp_value_to_double(v.get_internal());
  }

  Integer numerator(const Value& v) {
    Integer res;
    lp_value_get_num(v.get_internal(), res.get_internal());
    return res;
  }
  Integer denominator(const Value& v) {
    Integer res;
    lp_value_get_den(v.get_internal(), res.get_internal());
    return res;
  }

  Integer ceil(const Value& v) {
    Integer res;
    lp_value_ceiling(v.get_internal(), res.get_internal());
    return res;
  }
  Integer floor(const Value& v) {
    Integer res;
    lp_value_floor(v.get_internal(), res.get_internal());
    return res;
  }

  Value value_between(const lp_value_t* lhs, bool l_strict,
                      const lp_value_t* rhs, bool r_strict) {
    Value res;
    lp_value_get_value_between(lhs, l_strict ? 1 : 0, rhs, r_strict ? 1 : 0,
                               res.get_internal());
    return res;
  }

  Value value_between(const Value& lhs, bool l_strict, const Value& rhs,
                      bool r_strict) {
    return value_between(lhs.get_internal(), l_strict, rhs.get_internal(),
                         r_strict);
  }

  int approximate_size(const Value& lower, const Value& upper) {
    return lp_value_get_distance_size_approx(lower.get_internal(),
                                             upper.get_internal());
  }

}  // namespace poly
