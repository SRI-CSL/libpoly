#include "polyxx/integer.h"

#include <cassert>
#include <iostream>
#include <utility>

#include "polyxx/rational.h"

namespace poly {

  Integer::Integer() { lp_integer_construct(&mInt); }

  Integer::Integer(int i) : Integer(IntegerRing::Z, i) {}
  Integer::Integer(long i) : Integer(IntegerRing::Z, i) {}
  Integer::Integer(const IntegerRing& ir, long i) {
    lp_integer_construct_from_int(ir.get_internal(), &mInt, i);
  }

  Integer::Integer(const char* x, int base)
      : Integer(IntegerRing::Z, x, base) {}
  Integer::Integer(const IntegerRing& ir, const char* x, int base) {
    lp_integer_construct_from_string(ir.get_internal(), &mInt, x, base);
  }

  Integer::Integer(const Integer& i) : Integer(IntegerRing::Z, i) {}
  Integer::Integer(const IntegerRing& ir, const Integer& i) {
    lp_integer_construct_copy(ir.get_internal(), &mInt, i.get_internal());
  }

  Integer::Integer(const Rational& r) : Integer(IntegerRing::Z, r) {}
  Integer::Integer(const IntegerRing& ir, const Rational& r) {
    lp_integer_construct_from_rational(ir.get_internal(), &mInt,
                                       r.get_internal());
  }

  Integer::Integer(const mpz_class& m) : Integer(IntegerRing::Z, m) {}
  Integer::Integer(const IntegerRing& ir, const mpz_class& m) {
    lp_integer_construct_copy(ir.get_internal(), &mInt, m.get_mpz_t());
  }

  Integer::Integer(const lp_integer_t* i) : Integer(IntegerRing::Z, i) {}
  Integer::Integer(const IntegerRing& ir, const lp_integer_t* i) {
    lp_integer_construct_copy(ir.get_internal(), &mInt, i);
  }

  Integer::~Integer() { lp_integer_destruct(&mInt); }

  Integer& Integer::operator=(const Integer& i) {
    return assign(IntegerRing::Z, i);
  }
  Integer& Integer::assign(const IntegerRing& ir, const Integer& i) {
    lp_integer_assign(ir.get_internal(), &mInt, i.get_internal());
    return *this;
  }
  Integer& Integer::operator=(Integer&& i) {
    return assign(IntegerRing::Z, std::move(i));
  }
  Integer& Integer::assign(const IntegerRing& ir, Integer&& i) {
    lp_integer_assign(ir.get_internal(), &mInt, i.get_internal());
    return *this;
  }
  Integer& Integer::operator=(long i) { return assign(IntegerRing::Z, i); }
  Integer& Integer::assign(const IntegerRing& ir, long i) {
    lp_integer_assign_int(ir.get_internal(), &mInt, i);
    return *this;
  }

  lp_integer_t* Integer::get_internal() { return &mInt; }

  const lp_integer_t* Integer::get_internal() const { return &mInt; }

  std::ostream& operator<<(std::ostream& os, const Integer& i) {
    return stream_ptr(os, lp_integer_to_string(i.get_internal()));
  }

  std::size_t bit_size(const Integer& i) {
    return lp_integer_bits(i.get_internal());
  }

  bool operator==(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) == 0;
  }
  bool operator==(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) == 0;
  }
  bool operator==(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) == 0;
  }
  bool operator!=(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) != 0;
  }
  bool operator!=(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) != 0;
  }
  bool operator!=(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) != 0;
  }
  bool operator<(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) < 0;
  }
  bool operator<(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) < 0;
  }
  bool operator<(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) < 0;
  }
  bool operator<=(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) <= 0;
  }
  bool operator<=(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) <= 0;
  }
  bool operator<=(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) <= 0;
  }
  bool operator>(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) > 0;
  }
  bool operator>(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) > 0;
  }
  bool operator>(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) > 0;
  }
  bool operator>=(const Integer& lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) >= 0;
  }
  bool operator>=(const Integer& lhs, long rhs) {
    return compare(IntegerRing::Z, lhs, rhs) >= 0;
  }
  bool operator>=(long lhs, const Integer& rhs) {
    return compare(IntegerRing::Z, lhs, rhs) >= 0;
  }

  int compare(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    return lp_integer_cmp(ir.get_internal(), lhs.get_internal(),
                          rhs.get_internal());
  }
  int compare(const IntegerRing& ir, const Integer& lhs, long rhs) {
    return lp_integer_cmp_int(ir.get_internal(), lhs.get_internal(), rhs);
  }
  int compare(const IntegerRing& ir, long lhs, const Integer& rhs) {
    return -lp_integer_cmp_int(ir.get_internal(), rhs.get_internal(), lhs);
  }

  bool divides(const Integer& lhs, const Integer& rhs) {
    return divides(IntegerRing::Z, lhs, rhs);
  }
  bool divides(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    return lp_integer_divides(ir.get_internal(), lhs.get_internal(),
                              rhs.get_internal());
  }

  void swap(Integer& lhs, Integer& rhs) {
    lp_integer_swap(lhs.get_internal(), rhs.get_internal());
  }

  Integer& operator++(Integer& i) { return increment(IntegerRing::Z, i); }
  Integer& operator--(Integer& i) { return decrement(IntegerRing::Z, i); }
  Integer operator++(Integer& i, int) {
    Integer res(i);
    ++i;
    return res;
  }
  Integer operator--(Integer& i, int) {
    Integer res(i);
    --i;
    return res;
  }

  Integer& increment(const IntegerRing& ir, Integer& i) {
    lp_integer_inc(ir.get_internal(), i.get_internal());
    return i;
  }
  Integer& decrement(const IntegerRing& ir, Integer& i) {
    lp_integer_dec(ir.get_internal(), i.get_internal());
    return i;
  }

  Integer operator+(const Integer& lhs, const Integer& rhs) {
    return add(IntegerRing::Z, lhs, rhs);
  }
  Integer& operator+=(Integer& lhs, const Integer& rhs) {
    return add_assign(IntegerRing::Z, lhs, rhs);
  }
  Integer add(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    Integer res(lhs);
    return add_assign(ir, res, rhs);
  }
  Integer& add_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs) {
    lp_integer_add(ir.get_internal(), lhs.get_internal(), lhs.get_internal(),
                   rhs.get_internal());
    return lhs;
  }

  Integer operator-(const Integer& lhs, const Integer& rhs) {
    return sub(IntegerRing::Z, lhs, rhs);
  }
  Integer& operator-=(Integer& lhs, const Integer& rhs) {
    return sub_assign(IntegerRing::Z, lhs, rhs);
  }
  Integer sub(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    Integer res(lhs);
    return sub_assign(ir, res, rhs);
  }

  Integer& sub_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs) {
    lp_integer_sub(ir.get_internal(), lhs.get_internal(), lhs.get_internal(),
                   rhs.get_internal());
    return lhs;
  }

  Integer operator-(const Integer& i) { return neg(IntegerRing::Z, i); }
  Integer neg(const IntegerRing& ir, const Integer& i) {
    Integer res;
    lp_integer_neg(ir.get_internal(), res.get_internal(), i.get_internal());
    return res;
  }

  Integer abs(const Integer& i) { return abs(IntegerRing::Z, i); }
  Integer abs(const IntegerRing& ir, const Integer& i) {
    Integer res;
    lp_integer_abs(ir.get_internal(), res.get_internal(), i.get_internal());
    return res;
  }

  Integer inverse(const IntegerRing& ir, const Integer& i) {
    Integer res;
    lp_integer_inv(ir.get_internal(), res.get_internal(), i.get_internal());
    return res;
  }

  Integer operator*(const Integer& lhs, const Integer& rhs) {
    return mul(IntegerRing::Z, lhs, rhs);
  }
  Integer operator*(const Integer& lhs, long rhs) {
    return mul(IntegerRing::Z, lhs, rhs);
  }
  Integer operator*(long lhs, const Integer& rhs) {
    return mul(IntegerRing::Z, lhs, rhs);
  }
  Integer& operator*=(Integer& lhs, const Integer& rhs) {
    return mul_assign(IntegerRing::Z, lhs, rhs);
  }
  Integer& operator*=(Integer& lhs, long rhs) {
    return mul_assign(IntegerRing::Z, lhs, rhs);
  }

  Integer mul(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    Integer res(lhs);
    return mul_assign(ir, res, rhs);
  }
  Integer mul(const IntegerRing& ir, const Integer& lhs, long rhs) {
    Integer res(lhs);
    return mul_assign(ir, res, rhs);
  }
  Integer mul(const IntegerRing& ir, long lhs, const Integer& rhs) {
    return mul(ir, rhs, lhs);
  }
  Integer& mul_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs) {
    lp_integer_mul(ir.get_internal(), lhs.get_internal(), lhs.get_internal(),
                   rhs.get_internal());
    return lhs;
  }
  Integer& mul_assign(const IntegerRing& ir, Integer& lhs, long rhs) {
    lp_integer_mul_int(ir.get_internal(), lhs.get_internal(),
                       lhs.get_internal(), rhs);
    return lhs;
  }

  Integer mul_pow2(const Integer& lhs, unsigned rhs) {
    return mul_pow2(IntegerRing::Z, lhs, rhs);
  }
  Integer mul_pow2(const IntegerRing& ir, const Integer& lhs, unsigned rhs) {
    Integer res;
    lp_integer_mul_pow2(ir.get_internal(), res.get_internal(),
                        lhs.get_internal(), rhs);
    return res;
  }

  Integer pow(const Integer& lhs, unsigned rhs) {
    return pow(IntegerRing::Z, lhs, rhs);
  }
  Integer pow(const IntegerRing& ir, const Integer& lhs, unsigned rhs) {
    Integer res;
    lp_integer_pow(ir.get_internal(), res.get_internal(), lhs.get_internal(),
                   rhs);
    return res;
  }

  Integer sqrt(const Integer& i) {
    assert(i >= Integer());
    Integer res;
    lp_integer_sqrt_Z(res.get_internal(), i.get_internal());
    return res;
  }

  Integer& add_mul(Integer& lhs, const Integer& a, const Integer& b) {
    return add_mul(IntegerRing::Z, lhs, a, b);
  }
  Integer& add_mul(const IntegerRing& ir, Integer& lhs, const Integer& a,
                   const Integer& b) {
    lp_integer_add_mul(ir.get_internal(), lhs.get_internal(), a.get_internal(),
                       b.get_internal());
    return lhs;
  }
  Integer& add_mul(Integer& lhs, const Integer& a, int b) {
    return add_mul(IntegerRing::Z, lhs, a, b);
  }
  Integer& add_mul(const IntegerRing& ir, Integer& lhs, const Integer& a, int b) {
    lp_integer_add_mul_int(ir.get_internal(), lhs.get_internal(),
                           a.get_internal(), b);
    return lhs;
  }

  Integer& sub_mul(Integer& lhs, const Integer& a, const Integer& b) {
    return sub_mul(IntegerRing::Z, lhs, a, b);
  }
  Integer& sub_mul(const IntegerRing& ir, Integer& lhs, const Integer& a,
                   const Integer& b) {
    lp_integer_sub_mul(ir.get_internal(), lhs.get_internal(), a.get_internal(),
                       b.get_internal());
    return lhs;
  }

  Integer operator/(const Integer& lhs, const Integer& rhs) {
    Integer res(lhs);
    return res /= rhs;
  }
  Integer& operator/=(Integer& lhs, const Integer& rhs) {
    lp_integer_div_Z(lhs.get_internal(), lhs.get_internal(),
                     rhs.get_internal());
    return lhs;
  }

  Integer operator%(const Integer& lhs, const Integer& rhs) {
    Integer res(lhs);
    return res %= rhs;
  }
  Integer& operator%=(Integer& lhs, const Integer& rhs) {
    lp_integer_rem_Z(lhs.get_internal(), lhs.get_internal(),
                     rhs.get_internal());
    return lhs;
  }

  Integer div_exact(const Integer& lhs, const Integer& rhs) {
    return div_exact(IntegerRing::Z, lhs, rhs);
  }
  Integer div_exact(const IntegerRing& ir, const Integer& lhs, const Integer& rhs) {
    Integer res;
    lp_integer_div_exact(ir.get_internal(), res.get_internal(),
                         lhs.get_internal(), rhs.get_internal());
    return res;
  }

  Integer div_rem(Integer& rem, const Integer& lhs, const Integer& rhs) {
    Integer res;
    lp_integer_div_rem_Z(res.get_internal(), rem.get_internal(),
                         lhs.get_internal(), rhs.get_internal());
    return res;
  }

  Integer div_rem_pow2(Integer& rem, const Integer& lhs, unsigned rhs) {
    Integer res;
    lp_integer_div_rem_pow2_Z(res.get_internal(), rem.get_internal(),
                              lhs.get_internal(), rhs);
    return res;
  }

  long to_int(const Integer& i) { return lp_integer_to_int(i.get_internal()); }

  double to_double(const Integer& i) {
    return lp_integer_to_double(i.get_internal());
  }

  bool is_prime(const Integer& i) {
    return lp_integer_is_prime(i.get_internal()) > 0;
  }

  bool is_zero(const Integer& i) { return is_zero(IntegerRing::Z, i); }
  bool is_zero(const IntegerRing& ir, const Integer& i) {
    return lp_integer_is_zero(ir.get_internal(), i.get_internal());
  }

  bool is_in_ring(const IntegerRing& ir, const Integer& i) {
    return lp_integer_in_ring(ir.get_internal(), i.get_internal());
  }

  std::size_t hash(const Integer& i) {
    return lp_integer_hash(i.get_internal());
  }

  int sgn(const Integer& i) { return sgn(IntegerRing::Z, i); }
  int sgn(const IntegerRing& ir, const Integer& i) {
    return lp_integer_sgn(ir.get_internal(), i.get_internal());
  }

  Integer gcd(const Integer& a, const Integer& b) {
    Integer res;
    lp_integer_gcd_Z(res.get_internal(), a.get_internal(), b.get_internal());
    return res;
  }
  Integer lcm(const Integer& a, const Integer& b) {
    Integer res;
    lp_integer_lcm_Z(res.get_internal(), a.get_internal(), b.get_internal());
    return res;
  }

}  // namespace poly
