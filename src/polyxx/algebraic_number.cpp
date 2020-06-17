#include "polyxx/algebraic_number.h"

#include <iostream>

namespace poly {

  AlgebraicNumber::AlgebraicNumber(const lp_algebraic_number_t* an) {
    lp_algebraic_number_construct_copy(get_internal(), an);
  }
  AlgebraicNumber::AlgebraicNumber() {
    lp_algebraic_number_construct_zero(get_internal());
  }
  AlgebraicNumber::AlgebraicNumber(const AlgebraicNumber& an) {
    lp_algebraic_number_construct_copy(get_internal(), an.get_internal());
  }
  AlgebraicNumber::AlgebraicNumber(AlgebraicNumber&& an) : AlgebraicNumber() {
    swap(*this, an);
  }
  AlgebraicNumber::AlgebraicNumber(const DyadicRational& dr) {
    lp_algebraic_number_construct_from_dyadic_rational(get_internal(),
                                                       dr.get_internal());
  }

  AlgebraicNumber::AlgebraicNumber(UPolynomial&& poly,
                                   const DyadicInterval& di) {
    lp_algebraic_number_construct(get_internal(), poly.release(),
                                  di.get_internal());
  }
  AlgebraicNumber::AlgebraicNumber(const UPolynomial& poly,
                                   const DyadicInterval& di) {
    lp_algebraic_number_construct(get_internal(), UPolynomial(poly).release(),
                                  di.get_internal());
  }
  AlgebraicNumber::~AlgebraicNumber() {
    lp_algebraic_number_destruct(get_internal());
  }
  AlgebraicNumber& AlgebraicNumber::operator=(const AlgebraicNumber& an) {
    lp_algebraic_number_destruct(get_internal());
    lp_algebraic_number_construct_copy(get_internal(), an.get_internal());
    return *this;
  }
  AlgebraicNumber& AlgebraicNumber::operator=(AlgebraicNumber&& an) {
    swap(*this, an);
    return *this;
  }

  lp_algebraic_number_t* AlgebraicNumber::get_internal() { return &mValue; }
  const lp_algebraic_number_t* AlgebraicNumber::get_internal() const {
    return &mValue;
  }

  std::ostream& operator<<(std::ostream& os, const AlgebraicNumber& v) {
    return stream_ptr(os, lp_algebraic_number_to_string(v.get_internal()));
  }

  void swap(AlgebraicNumber& lhs, AlgebraicNumber& rhs) {
    lp_algebraic_number_swap(lhs.get_internal(), rhs.get_internal());
  }
  int sgn(const AlgebraicNumber& an) {
    return lp_algebraic_number_sgn(an.get_internal());
  }

  bool operator==(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const AlgebraicNumber& lhs, const AlgebraicNumber& rhs) {
    return lp_algebraic_number_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) == 0;
  }
  bool operator!=(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) != 0;
  }
  bool operator<(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) < 0;
  }
  bool operator<=(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) <= 0;
  }
  bool operator>(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) > 0;
  }
  bool operator>=(const AlgebraicNumber& lhs, const Integer& rhs) {
    return lp_algebraic_number_cmp_integer(lhs.get_internal(),
                                           rhs.get_internal()) >= 0;
  }

  bool operator==(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Integer& lhs, const AlgebraicNumber& rhs) {
    return rhs <= lhs;
  }

  bool operator==(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) == 0;
  }
  bool operator!=(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) != 0;
  }
  bool operator<(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) < 0;
  }
  bool operator<=(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) <= 0;
  }
  bool operator>(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) > 0;
  }
  bool operator>=(const AlgebraicNumber& lhs, const DyadicRational& rhs) {
    return lp_algebraic_number_cmp_dyadic_rational(lhs.get_internal(),
                                                   rhs.get_internal()) >= 0;
  }

  bool operator==(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs != lhs;
  }
  bool operator<(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const DyadicRational& lhs, const AlgebraicNumber& rhs) {
    return rhs <= lhs;
  }

  bool operator==(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) == 0;
  }
  bool operator!=(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) != 0;
  }
  bool operator<(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) < 0;
  }
  bool operator<=(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) <= 0;
  }
  bool operator>(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) > 0;
  }
  bool operator>=(const AlgebraicNumber& lhs, const Rational& rhs) {
    return lp_algebraic_number_cmp_rational(lhs.get_internal(),
                                            rhs.get_internal()) >= 0;
  }

  bool operator==(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Rational& lhs, const AlgebraicNumber& rhs) {
    return rhs <= lhs;
  }

  double to_double(const AlgebraicNumber& an) {
    return lp_algebraic_number_to_double(an.get_internal());
  }

  Rational to_rational_approximation(const AlgebraicNumber& an) {
    Rational res;
    lp_algebraic_number_to_rational(an.get_internal(), res.get_internal());
    return res;
  }

  const DyadicRational& get_lower_bound(const AlgebraicNumber& an) {
    return get_lower(*detail::cast_from(&an.get_internal()->I));
  }
  const DyadicRational& get_upper_bound(const AlgebraicNumber& an) {
    return get_upper(*detail::cast_from(&an.get_internal()->I));
  }

  DyadicRational midpoint_dyadic(const AlgebraicNumber& an) {
    DyadicRational res;
    lp_algebraic_number_get_dyadic_midpoint(an.get_internal(),
                                            res.get_internal());
    return res;
  }
  Rational midpoint_rational(const AlgebraicNumber& an) {
    Rational res;
    lp_algebraic_number_get_rational_midpoint(an.get_internal(),
                                              res.get_internal());
    return res;
  }
  void refine(AlgebraicNumber& an) {
    lp_algebraic_number_refine(an.get_internal());
  }
  void refine_const(const AlgebraicNumber& an) {
    lp_algebraic_number_refine_const(an.get_internal());
  }

  AlgebraicNumber operator+(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs) {
    AlgebraicNumber res;
    lp_algebraic_number_add(res.get_internal(), lhs.get_internal(),
                            rhs.get_internal());
    return res;
  }
  AlgebraicNumber operator-(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs) {
    AlgebraicNumber res;
    lp_algebraic_number_sub(res.get_internal(), lhs.get_internal(),
                            rhs.get_internal());
    return res;
  }

  AlgebraicNumber operator-(const AlgebraicNumber& an) {
    AlgebraicNumber res;
    lp_algebraic_number_neg(res.get_internal(), an.get_internal());
    return res;
  }
  AlgebraicNumber operator*(const AlgebraicNumber& lhs,
                            const AlgebraicNumber& rhs) {
    AlgebraicNumber res;
    lp_algebraic_number_mul(res.get_internal(), lhs.get_internal(),
                            rhs.get_internal());
    return res;
  }
  AlgebraicNumber pow(const AlgebraicNumber& lhs, unsigned n) {
    AlgebraicNumber res;
    lp_algebraic_number_pow(res.get_internal(), lhs.get_internal(), n);
    return res;
  }

  bool is_rational(const AlgebraicNumber& an) {
    return lp_algebraic_number_is_rational(an.get_internal());
  }
  bool is_integer(const AlgebraicNumber& an) {
    return lp_algebraic_number_is_integer(an.get_internal());
  }
  bool is_zero(const AlgebraicNumber& an) { return an == AlgebraicNumber(); }
  bool is_one(const AlgebraicNumber& an) {
    return an == AlgebraicNumber(DyadicRational(1));
  }

  Integer ceil(const AlgebraicNumber& an) {
    Integer res;
    lp_algebraic_number_ceiling(an.get_internal(), res.get_internal());
    return res;
  }
  Integer floor(const AlgebraicNumber& an) {
    Integer res;
    lp_algebraic_number_floor(an.get_internal(), res.get_internal());
    return res;
  }

}  // namespace poly
