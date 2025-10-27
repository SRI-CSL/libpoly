#pragma once

#include <gmpxx.h>

#include <iosfwd>

#include "../integer.h"
#include "integer_ring.h"

namespace poly {

  class Rational;

  /**
   * Implements a wrapper for lp_integer_t.
   */
  class Integer {
    /** The actual integer. */
    lp_integer_t mInt;

   public:
    /** Constructs zero. */
    Integer();
    /** Constructs from an int. */
    Integer(int i);
    /** Constructs from a long. */
    Integer(long i);
    /** Constructs from a long into the given ring. */
    Integer(const IntegerRing& ir, long i);
    /** Constructs from a string. */
    Integer(const char* x, int base);
    /** Constructs from a string into the given ring. */
    Integer(const IntegerRing& ir, const char* x, int base);
    /** Constructs as copy. */
    Integer(const Integer& i);
    /** Constructs as copy into the given ring. */
    Integer(const IntegerRing& ir, const Integer& i);
    /** Constructs as copy, assuming the rational is indeed an integer. */
    explicit Integer(const Rational& r);
    /** Constructs as copy, assuming the rational is indeed an integer, into the
     * given ring. */
    Integer(const IntegerRing& ir, const Rational& r);

    /** Construct from a mpz_class, which is the underlying representation
     * anyway. */
    explicit Integer(const mpz_class& m);
    /** Construct from a mpz_class, which is the underlying representation
     * anyway, into the given ring. */
    Integer(const IntegerRing& ir, const mpz_class& m);
    /** Construct from an internal lp_integer_t pointer. */
    explicit Integer(const lp_integer_t* i);
    /** Construct from an internal lp_integer_t pointer into the given ring. */
    Integer(const IntegerRing& ir, const lp_integer_t* i);

    /** Custom destructor. */
    ~Integer();
    /** Assign from an Integer. */
    Integer& operator=(const Integer& i);
    /** Assign from an Integer into the given ring. */
    Integer& assign(const IntegerRing& ir, const Integer& i);
    /** Move from an Integer. */
    Integer& operator=(Integer&& i);
    /** Move from an Integer into the given ring. */
    Integer& assign(const IntegerRing& ir, Integer&& i);
    /** Assign from the given integer. */
    Integer& operator=(long i);
    /** Assign from the given integer into the given ring. */
    Integer& assign(const IntegerRing& ir, long i);

    /** Get a non-const pointer to the internal lp_integer_t. Handle with care!
     */
    lp_integer_t* get_internal();
    /** Get a const pointer to the internal lp_integer_t. */
    const lp_integer_t* get_internal() const;
  };

  /** Make sure that we can cast between Integer and lp_integer_t. */
  static_assert(sizeof(Integer) == sizeof(lp_integer_t),
                "Please check the size of Integer.");
  static_assert(sizeof(Integer) == sizeof(mpz_class),
                "Please check the size of Integer.");
  namespace detail {
    /** Non-const cast from an Integer to a lp_integer_t. */
    inline lp_integer_t* cast_to(Integer* i) {
      return reinterpret_cast<lp_integer_t*>(i);
    }
    /** Const cast from an Integer to a lp_integer_t. */
    inline const lp_integer_t* cast_to(const Integer* i) {
      return reinterpret_cast<const lp_integer_t*>(i);
    }
    /** Non-const cast from an Integer to a mpz_class. */
    inline mpz_class* cast_to_gmp(Integer* i) {
      return reinterpret_cast<mpz_class*>(i);
    }
    /** Const cast from an Integer to a mpz_class. */
    inline const mpz_class* cast_to_gmp(const Integer* i) {
      return reinterpret_cast<const mpz_class*>(i);
    }
    /** Non-const cast from a lp_integer_t to an Integer. */
    inline Integer* cast_from(lp_integer_t* i) {
      return reinterpret_cast<Integer*>(i);
    }
    /** Const cast from a lp_integer_t to an Integer. */
    inline const Integer* cast_from(const lp_integer_t* i) {
      return reinterpret_cast<const Integer*>(i);
    }
    /** Non-const cast from a mpz_class to an Integer. */
    inline Integer* cast_from(mpz_class* i) {
      return reinterpret_cast<Integer*>(i);
    }
    /** Const cast from a mpz_class to an Integer. */
    inline const Integer* cast_from(const mpz_class* i) {
      return reinterpret_cast<const Integer*>(i);
    }
  }  // namespace detail

  /** Stream the given Integer to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Integer& i);

  /** Number of bits needed for the given integer. */
  std::size_t bit_size(const Integer& i);

  /** Compare two integers. */
  bool operator==(const Integer& lhs, const Integer& rhs);
  bool operator==(const Integer& lhs, long rhs);
  bool operator==(long lhs, const Integer& rhs);
  /** Compare two integers. */
  bool operator!=(const Integer& lhs, const Integer& rhs);
  bool operator!=(const Integer& lhs, long rhs);
  bool operator!=(long lhs, const Integer& rhs);
  /** Compare two integers according to the lexicographic ordering on (lower bound,upper bound). */
  bool operator<(const Integer& lhs, const Integer& rhs);
  bool operator<(const Integer& lhs, long rhs);
  bool operator<(long lhs, const Integer& rhs);
  /** Compare two integers according to the lexicographic ordering on (lower bound,upper bound). */
  bool operator<=(const Integer& lhs, const Integer& rhs);
  bool operator<=(const Integer& lhs, long rhs);
  bool operator<=(long lhs, const Integer& rhs);
  /** Compare two integers according to the lexicographic ordering on (lower bound,upper bound). */
  bool operator>(const Integer& lhs, const Integer& rhs);
  bool operator>(const Integer& lhs, long rhs);
  bool operator>(long lhs, const Integer& rhs);
  /** Compare two integers according to the lexicographic ordering on (lower bound,upper bound). */
  bool operator>=(const Integer& lhs, const Integer& rhs);
  bool operator>=(const Integer& lhs, long rhs);
  bool operator>=(long lhs, const Integer& rhs);

  /** Compare two integers over the given ring. */
  int compare(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);
  /** Compare two integers over the given ring. */
  int compare(const IntegerRing& ir, const Integer& lhs, long rhs);
  /** Compare two integers over the given ring. */
  int compare(const IntegerRing& ir, long lhs, const Integer& rhs);

  /** Checks whether lhs divides rhs. */
  bool divides(const Integer& lhs, const Integer& rhs);
  /** Checks whether lhs divides rhs over the given ring. */
  bool divides(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);

  /** Swaps the contents of two integers. */
  void swap(Integer& lhs, Integer& rhs);

  /** Pre-increment an integer. */
  Integer& operator++(Integer& i);
  /** Pre-decrement an integer. */
  Integer& operator--(Integer& i);
  /** Post-increment an integer. */
  Integer operator++(Integer& i, int);
  /** Post-decrement an integer. */
  Integer operator--(Integer& i, int);
  /** Pre-increment an integer. */
  Integer& increment(const IntegerRing& ir, Integer& i);
  /** Pre-decrement an integer. */
  Integer& decrement(const IntegerRing& ir, Integer& i);

  /** Add two integers. */
  Integer operator+(const Integer& lhs, const Integer& rhs);
  /** Add and assign two integers. */
  Integer& operator+=(Integer& lhs, const Integer& rhs);
  /** Add two integers in the given ring. */
  Integer add(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);
  /** Add and assign two integers in the given ring. */
  Integer& add_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs);

  /** Subtract two integers. */
  Integer operator-(const Integer& lhs, const Integer& rhs);
  /** Subtract and assign two integers. */
  Integer& operator-=(Integer& lhs, const Integer& rhs);
  /** Subtract two integers in the given ring. */
  Integer sub(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);
  /** Subtract and assign two integers in the given ring. */
  Integer& sub_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs);

  /** Negate an integer. */
  Integer operator-(const Integer& i);
  /** Negate an integer in the given ring. */
  Integer neg(const IntegerRing& ir, const Integer& i);

  /** Compute the absolute value. */
  Integer abs(const Integer& i);
  /** Compute the absolute value in the given ring. */
  Integer abs(const IntegerRing& ir, const Integer& i);

  /** Compute the inverse in the given ring.
   * Note that inverses do not exist in Z (except for -1 and 1).
   */
  Integer inverse(const IntegerRing& ir, const Integer& i);

  /** Multiply two integers. */
  Integer operator*(const Integer& lhs, const Integer& rhs);
  /** Multiply two integers. */
  Integer operator*(const Integer& lhs, long rhs);
  /** Multiply two integers. */
  Integer operator*(long lhs, const Integer& rhs);
  /** Multiply and assign two integers. */
  Integer& operator*=(Integer& lhs, const Integer& rhs);
  /** Multiply and assign two integers. */
  Integer& operator*=(Integer& lhs, long rhs);
  /** Multiply two integers in the given ring. */
  Integer mul(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);
  /** Multiply two integers in the given ring. */
  Integer mul(const IntegerRing& ir, const Integer& lhs, long rhs);
  /** Multiply two integers in the given ring. */
  Integer mul(const IntegerRing& ir, long lhs, const Integer& rhs);
  /** Multiply and assign two integers in the given ring. */
  Integer& mul_assign(const IntegerRing& ir, Integer& lhs, const Integer& rhs);
  /** Multiply and assign two integers in the given ring. */
  Integer& mul_assign(const IntegerRing& ir, Integer& lhs, long rhs);

  /** Compute lhs * 2^rhs. */
  Integer mul_pow2(const Integer& lhs, unsigned rhs);
  /** Compute lhs * 2^rhs in the given ring. */
  Integer mul_pow2(const IntegerRing& ir, const Integer& lhs, unsigned rhs);

  /** Compute lhs^rhs. */
  Integer pow(const Integer& lhs, unsigned rhs);
  /** Compute lhs^rhs in the given ring. */
  Integer pow(const IntegerRing& ir, const Integer& lhs, unsigned rhs);

  /** Compute the (truncated part of the) square root. */
  Integer sqrt(const Integer& i);

  /** Compute lhs += a * b. */
  Integer& add_mul(Integer& lhs, const Integer& a, const Integer& b);
  /** Compute lhs += a * b in the given ring. */
  Integer& add_mul(const IntegerRing& ir, Integer& lhs, const Integer& a,
                   const Integer& b);
  /** Compute lhs += a * b. */
  Integer& add_mul(Integer& lhs, const Integer& a, int b);
  /** Compute lhs += a * b in the given ring. */
  Integer& add_mul(const IntegerRing& ir, Integer& lhs, const Integer& a, int b);

  /** Compute lhs -= a * b. */
  Integer& sub_mul(Integer& lhs, const Integer& a, const Integer& b);
  /** Compute lhs -= a * b in the given ring. */
  Integer& sub_mul(const IntegerRing& ir, Integer& lhs, const Integer& a,
                   const Integer& b);

  /** Compute the (truncated part of the) quotient of two integers. */
  Integer operator/(const Integer& lhs, const Integer& rhs);
  /** Compute and assign the (truncated part of the) quotient of two integers.
   */
  Integer& operator/=(Integer& lhs, const Integer& rhs);

  /** Compute the remainder of two integers. */
  Integer operator%(const Integer& lhs, const Integer& rhs);
  /** Compute and assign the remainder of two integers. */
  Integer& operator%=(Integer& lhs, const Integer& rhs);

  /** Compute the quotient of two integers, assuming that divides(lhs, rhs). */
  Integer div_exact(const Integer& lhs, const Integer& rhs);
  /** Compute the quotient of two integers in the given ring, assuming that
   * divides(lhs, rhs). */
  Integer div_exact(const IntegerRing& ir, const Integer& lhs, const Integer& rhs);

  /** Compute the quotient and remainder of two integers. */
  Integer div_rem(Integer& rem, const Integer& lhs, const Integer& rhs);

  /** Compute the quotient and remainder of lhs and 2^rhs. */
  Integer div_rem_pow2(Integer& rem, const Integer& lhs, unsigned rhs);

  /** Returns the integer as a long. May get truncated, but keeps the correct
   * sign. */
  long to_int(const Integer& i);
  /** Returns the integer as a double. */
  double to_double(const Integer& i);

  /** Checks whether the integer is a prime number.
   * If the integer is negative, checks whether abs(i) is a prime number.
   */
  bool is_prime(const Integer& i);
  /** Checks whether the integer is zero. */
  bool is_zero(const Integer& i);
  /** Checks whether the integer is zero in the given ring. */
  bool is_zero(const IntegerRing& ir, const Integer& i);
  /** Checks whether the integer is in the given ring. */
  bool is_in_ring(const IntegerRing& ir, const Integer& i);

  /** Computes the hash of an integer. */
  std::size_t hash(const Integer& i);

  /** Computes the sign of an integer. */
  int sgn(const Integer& i);
  /** Computes the sign of an integer in the given ring. */
  int sgn(const IntegerRing& ir, const Integer& i);

  /** Computes the greatest common divisor of two integers. */
  Integer gcd(const Integer& a, const Integer& b);
  /** Computes the least common multiple of two integers. */
  Integer lcm(const Integer& a, const Integer& b);

}  // namespace poly

namespace std {
  template <>
  struct hash<poly::Integer> {
    std::size_t operator()(const poly::Integer& i) const {
      return poly::hash(i);
    }
  };
}  // namespace std
