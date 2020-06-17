#include "polyxx/rational.h"

#include <iostream>

namespace poly {

  Rational::Rational() { lp_rational_construct(&mRat); }
  Rational::Rational(int i) { lp_rational_construct_from_int(&mRat, i, 1); }

  Rational::Rational(const Rational& r) {
    lp_rational_construct_copy(&mRat, r.get_internal());
  }
  Rational::Rational(Rational&& r) {
    lp_rational_construct_copy(&mRat, r.get_internal());
  }

  Rational::Rational(const Integer& num, const Integer& denom) {
    lp_rational_construct_from_div(&mRat, num.get_internal(),
                                   denom.get_internal());
  }
  Rational::Rational(long num, unsigned long denom) {
    lp_rational_construct_from_int(&mRat, num, denom);
  }
  Rational::Rational(const Integer& i) {
    lp_rational_construct_from_integer(&mRat, i.get_internal());
  }
  Rational::Rational(double d) { lp_rational_construct_from_double(&mRat, d); }

  Rational::Rational(const mpq_class& m) {
    lp_rational_construct_copy(get_internal(), m.get_mpq_t());
  }

  Rational::~Rational() { lp_rational_destruct(&mRat); }

  Rational& Rational::operator=(const Rational& r) {
    lp_rational_assign(&mRat, r.get_internal());
    return *this;
  }
  Rational& Rational::operator=(Rational&& r) {
    lp_rational_assign(&mRat, r.get_internal());
    return *this;
  }
  lp_rational_t* Rational::get_internal() { return &mRat; }
  const lp_rational_t* Rational::get_internal() const { return &mRat; }

  std::ostream& operator<<(std::ostream& os, const Rational& r) {
    return stream_ptr(os, lp_rational_to_string(r.get_internal()));
  }

  double to_double(const Rational& r) {
    return lp_rational_to_double(r.get_internal());
  }

  int sgn(const Rational& r) { return lp_rational_sgn(r.get_internal()); }

  bool operator==(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const Rational& lhs, const Rational& rhs) {
    return lp_rational_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const Rational& lhs, const Integer& rhs) {
    return lp_rational_cmp_integer(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const Integer& lhs, const Rational& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Integer& lhs, const Rational& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Integer& lhs, const Rational& rhs) { return rhs > lhs; }
  bool operator<=(const Integer& lhs, const Rational& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Integer& lhs, const Rational& rhs) { return rhs < lhs; }
  bool operator>=(const Integer& lhs, const Rational& rhs) {
    return rhs <= lhs;
  }

  void swap(Rational& lhs, Rational& rhs) {
    lp_rational_swap(lhs.get_internal(), rhs.get_internal());
  }

  Rational operator+(const Rational& lhs, const Rational& rhs) {
    Rational res;
    lp_rational_add(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }
  Rational operator+(const Rational& lhs, const Integer& rhs) {
    Rational res;
    lp_rational_add_integer(res.get_internal(), lhs.get_internal(),
                            rhs.get_internal());
    return res;
  }
  Rational operator+(const Integer& lhs, const Rational& rhs) {
    return rhs + lhs;
  }

  Rational operator-(const Rational& lhs, const Rational& rhs) {
    Rational res;
    lp_rational_sub(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }
  Rational operator-(const Rational& r) {
    Rational res;
    lp_rational_neg(res.get_internal(), r.get_internal());
    return res;
  }
  Rational inverse(const Rational& r) {
    Rational res;
    lp_rational_inv(res.get_internal(), r.get_internal());
    return res;
  }

  Rational operator*(const Rational& lhs, const Rational& rhs) {
    Rational res;
    lp_rational_mul(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }

  Rational mul_2exp(const Rational& lhs, unsigned n) {
    Rational res;
    lp_rational_mul_2exp(res.get_internal(), lhs.get_internal(), n);
    return res;
  }

  Rational pow(const Rational& r, unsigned n) {
    Rational res;
    lp_rational_pow(res.get_internal(), r.get_internal(), n);
    return res;
  }

  Rational operator/(const Rational& lhs, const Rational& rhs) {
    Rational res;
    lp_rational_div(res.get_internal(), lhs.get_internal(), rhs.get_internal());
    return res;
  }

  Rational div_2exp(const Rational& lhs, unsigned n) {
    Rational res;
    lp_rational_div_2exp(res.get_internal(), lhs.get_internal(), n);
    return res;
  }

  const Integer& numerator(const Rational& r) {
    return *detail::cast_from(mpq_numref(r.get_internal()));
  }
  const Integer& denominator(const Rational& r) {
    return *detail::cast_from(mpq_denref(r.get_internal()));
  }

  bool is_integer(const Rational& r) {
    return lp_rational_is_integer(r.get_internal());
  }

  Integer ceil(const Rational& r) {
    Integer res;
    lp_rational_ceiling(r.get_internal(), res.get_internal());
    return res;
  }
  Integer floor(const Rational& r) {
    Integer res;
    lp_rational_floor(r.get_internal(), res.get_internal());
    return res;
  }

}  // namespace poly
