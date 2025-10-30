#include "polyxx/dyadic_rational.h"

#include <iostream>

#include "polyxx/rational.h"

namespace poly {

  DyadicRational::DyadicRational() { lp_dyadic_rational_construct(&mDRat); }

  DyadicRational::DyadicRational(const DyadicRational& r) {
    lp_dyadic_rational_construct_copy(&mDRat, r.get_internal());
  }
  DyadicRational::DyadicRational(DyadicRational&& r) {
    lp_dyadic_rational_construct_copy(&mDRat, r.get_internal());
  }

  DyadicRational::DyadicRational(long a, unsigned long n) {
    lp_dyadic_rational_construct_from_int(&mDRat, a, n);
  }
  DyadicRational::DyadicRational(const Integer& i) {
    lp_dyadic_rational_construct_from_integer(&mDRat, i.get_internal());
  }
  DyadicRational::DyadicRational(double d) {
    lp_dyadic_rational_construct_from_double(&mDRat, d);
  }
  DyadicRational::DyadicRational(int i) : DyadicRational(i, 0) {}
  DyadicRational::DyadicRational(long i) : DyadicRational(i, 0) {}

  DyadicRational::DyadicRational(const lp_dyadic_rational_t* dr) {
    lp_dyadic_rational_construct_copy(&mDRat, dr);
  }

  DyadicRational::~DyadicRational() { lp_dyadic_rational_destruct(&mDRat); }

  DyadicRational& DyadicRational::operator=(const DyadicRational& r) {
    lp_dyadic_rational_assign(&mDRat, r.get_internal());
    return *this;
  }
  DyadicRational& DyadicRational::operator=(DyadicRational&& r) {
    lp_dyadic_rational_assign(&mDRat, r.get_internal());
    return *this;
  }

  DyadicRational::operator Rational() const {
    Rational res;
    lp_rational_construct_from_dyadic(res.get_internal(), get_internal());
    return res;
  }

  lp_dyadic_rational_t* DyadicRational::get_internal() { return &mDRat; }
  const lp_dyadic_rational_t* DyadicRational::get_internal() const {
    return &mDRat;
  }

  std::ostream& operator<<(std::ostream& os, const DyadicRational& r) {
    return stream_ptr(os, lp_dyadic_rational_to_string(r.get_internal()));
  }

  double to_double(const DyadicRational& r) {
    return lp_dyadic_rational_to_double(r.get_internal());
  }

  int sgn(const DyadicRational& r) {
    return lp_dyadic_rational_sgn(r.get_internal());
  }

  bool operator==(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicRational& lhs, const DyadicRational& rhs) {
    return lp_dyadic_rational_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
  }

  bool operator==(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicRational& lhs, const Integer& rhs) {
    return lp_dyadic_rational_cmp_integer(lhs.get_internal(),
                                          rhs.get_internal()) >= 0;
  }

  bool operator==(const Integer& lhs, const DyadicRational& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Integer& lhs, const DyadicRational& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Integer& lhs, const DyadicRational& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Integer& lhs, const DyadicRational& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Integer& lhs, const DyadicRational& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Integer& lhs, const DyadicRational& rhs) {
    return rhs <= lhs;
  }

  bool operator==(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) == 0;
  }
  bool operator!=(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) != 0;
  }
  bool operator<(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) < 0;
  }
  bool operator<=(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) <= 0;
  }
  bool operator>(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) > 0;
  }
  bool operator>=(const DyadicRational& lhs, const Rational& rhs) {
    return lp_dyadic_rational_cmp_rational(lhs.get_internal(),
                                           rhs.get_internal()) >= 0;
  }

  bool operator==(const Rational& lhs, const DyadicRational& rhs) {
    return rhs == lhs;
  }
  bool operator!=(const Rational& lhs, const DyadicRational& rhs) {
    return rhs != lhs;
  }
  bool operator<(const Rational& lhs, const DyadicRational& rhs) {
    return rhs > lhs;
  }
  bool operator<=(const Rational& lhs, const DyadicRational& rhs) {
    return rhs >= lhs;
  }
  bool operator>(const Rational& lhs, const DyadicRational& rhs) {
    return rhs < lhs;
  }
  bool operator>=(const Rational& lhs, const DyadicRational& rhs) {
    return rhs <= lhs;
  }

  void swap(DyadicRational& lhs, DyadicRational& rhs) {
    lp_dyadic_rational_swap(lhs.get_internal(), rhs.get_internal());
  }

  DyadicRational operator+(const DyadicRational& lhs,
                           const DyadicRational& rhs) {
    DyadicRational res;
    lp_dyadic_rational_add(res.get_internal(), lhs.get_internal(),
                           rhs.get_internal());
    return res;
  }
  DyadicRational operator+(const DyadicRational& lhs, const Integer& rhs) {
    DyadicRational res;
    lp_dyadic_rational_add_integer(res.get_internal(), lhs.get_internal(),
                                   rhs.get_internal());
    return res;
  }
  DyadicRational operator+(const Integer& lhs, const DyadicRational& rhs) {
    return rhs + lhs;
  }

  DyadicRational operator-(const DyadicRational& lhs,
                           const DyadicRational& rhs) {
    DyadicRational res;
    lp_dyadic_rational_sub(res.get_internal(), lhs.get_internal(),
                           rhs.get_internal());
    return res;
  }
  DyadicRational operator-(const DyadicRational& r) {
    DyadicRational res;
    lp_dyadic_rational_neg(res.get_internal(), r.get_internal());
    return res;
  }

  DyadicRational operator*(const DyadicRational& lhs,
                           const DyadicRational& rhs) {
    DyadicRational res;
    lp_dyadic_rational_mul(res.get_internal(), lhs.get_internal(),
                           rhs.get_internal());
    return res;
  }

  DyadicRational mul_2exp(const DyadicRational& lhs, unsigned long n) {
    DyadicRational res;
    lp_dyadic_rational_mul_2exp(res.get_internal(), lhs.get_internal(), n);
    return res;
  }

  DyadicRational pow(const DyadicRational& r, unsigned long n) {
    DyadicRational res;
    lp_dyadic_rational_pow(res.get_internal(), r.get_internal(), n);
    return res;
  }

  DyadicRational div_2exp(const DyadicRational& lhs, unsigned long n) {
    DyadicRational res;
    lp_dyadic_rational_div_2exp(res.get_internal(), lhs.get_internal(), n);
    return res;
  }

  Integer numerator(const DyadicRational& r) {
    Integer res;
    lp_dyadic_rational_get_num(r.get_internal(), res.get_internal());
    return res;
  }
  Integer denominator(const DyadicRational& r) {
    Integer res;
    lp_dyadic_rational_get_den(r.get_internal(), res.get_internal());
    return res;
  }

  bool is_integer(const DyadicRational& r) {
    return lp_dyadic_rational_is_integer(r.get_internal());
  }

  Integer ceil(const DyadicRational& r) {
    Integer res;
    lp_dyadic_rational_ceiling(r.get_internal(), res.get_internal());
    return res;
  }
  Integer floor(const DyadicRational& r) {
    Integer res;
    lp_dyadic_rational_floor(r.get_internal(), res.get_internal());
    return res;
  }

}  // namespace poly
