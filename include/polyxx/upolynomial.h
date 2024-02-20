#pragma once

#include <iosfwd>
#include <utility>
#include <vector>

#include "../integer.h"
#include "../upolynomial.h"
#include "dyadic_rational.h"
#include "integer.h"
#include "integer_ring.h"
#include "rational.h"
#include "variable.h"

namespace poly {

  /**
   * Implements a wrapper for lp_upolynomial_t.
   */
  class UPolynomial {
    /** The actual univariate polynomial. */
    deleting_unique_ptr<lp_upolynomial_t> mPoly;

   public:
    /** Create from a lp_upolynomial_t pointer, claiming it's ownership. */
    explicit UPolynomial(lp_upolynomial_t* poly);
    /** Copy from an internal lp_upolynomial_t pointer. */
    explicit UPolynomial(const lp_upolynomial_t* poly);

    /** Create the zero polynomial. */
    explicit UPolynomial();
    /** Create a constant polynomial. */
    explicit UPolynomial(const Integer& i);
    /** Create a constant polynomial. */
    explicit UPolynomial(long i);

    /** Create from integer coefficients. */
    explicit UPolynomial(const std::vector<Integer>& coefficients);
    /** Create from integer coefficients into the given integer ring. */
    UPolynomial(const IntegerRing& ir,
                const std::vector<Integer>& coefficients);
    /** Create from integer coefficients. */
    explicit UPolynomial(const std::vector<long>& coefficients);
    /** Create from integer coefficients into the given integer ring. */
    UPolynomial(const IntegerRing& ir, const std::vector<long>& coefficients);

    /** Create from integer coefficients. */
    explicit UPolynomial(std::initializer_list<long> coefficients);
    /** Create from integer coefficients into the given integer ring. */
    UPolynomial(const IntegerRing& ir,
                std::initializer_list<long> coefficients);

    /** Construct c * x^degree. */
    UPolynomial(std::size_t degree, long c);
    /** Construct c * x^degree into the given integer ring. */
    UPolynomial(const IntegerRing& ir, std::size_t degree, long c);

    /** Copy from a polynomial. */
    UPolynomial(const UPolynomial& poly);
    /** Move from a polynomial. */
    UPolynomial(UPolynomial&& poly);
    /** Copy from a polynomial into the given integer ring. */
    UPolynomial(const IntegerRing& ir, const UPolynomial& poly);

    /** Copy from a polynomial. */
    UPolynomial& operator=(const UPolynomial& poly);
    /** Move from a polynomial. */
    UPolynomial& operator=(UPolynomial&& poly);
    /** Assign from and take ownership of an internal lp_upolynomial_t pointer.
     */
    UPolynomial& operator=(lp_upolynomial_t* poly);

    /** Get a non-const pointer to the internal lp_upolynomial_t. Handle with
     * care! */
    lp_upolynomial_t* get_internal();
    /** Get a const pointer to the internal lp_upolynomial_t. */
    const lp_upolynomial_t* get_internal() const;
    /** Release the lp_upolynomial_t pointer. This yields ownership of the
     * returned pointer. */
    lp_upolynomial_t* release();
  };

  /** Return the degree of a polynomial. */
  std::size_t degree(const UPolynomial& p);

  /** Get the leading coefficient. */
  const Integer& leading_coefficient(const UPolynomial& p);
  /** Get the constant coefficient. */
  const Integer& constant_coefficient(const UPolynomial& p);

  /** Get all coefficients. */
  std::vector<Integer> coefficients(const UPolynomial& p);

  /** Stream the given UPolynomial to an output stream. */
  std::ostream& operator<<(std::ostream& os, const UPolynomial& p);

  /** Check if the polynomial is zero. */
  bool is_zero(const UPolynomial& p);
  /** Check if the polynomial is one. */
  bool is_one(const UPolynomial& p);
  /** Check if the polynomial is monic.
   * A polynomial is monic if its leading coefficient is one.
   */
  bool is_monic(const UPolynomial& p);
  /** Check if the polynomial is primitive.
   * A polynomial is primitive if the gcd of all coefficients is one and the
   * leading coefficient is positive.
   */
  bool is_primitive(const UPolynomial& p);

  /** Evaluate a polynomial at an integer. */
  Integer evaluate_at(const UPolynomial& p, const Integer& i);
  /** Evaluate a polynomial at a rational. */
  Rational evaluate_at(const UPolynomial& p, const Rational& r);
  /** Evaluate a polynomial at a dyadic rational. */
  DyadicRational evaluate_at(const UPolynomial& p, const DyadicRational& dr);

  /** Determine the sign of a polynomial at an integer. */
  int sign_at(const UPolynomial& p, const Integer& i);
  /** Determine the sign of a polynomial at a rational. */
  int sign_at(const UPolynomial& p, const Rational& r);
  /** Determine the sign of a polynomial at a dyadic rational. */
  int sign_at(const UPolynomial& p, const DyadicRational& dr);

  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator==(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator!=(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator<(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator<=(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator>(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compares two polynomials (using a lexicographic order on the
   * coefficients). */
  bool operator>=(const UPolynomial& lhs, const UPolynomial& rhs);

  /** Compute p(-x). */
  UPolynomial subst_x_neg(const UPolynomial& p);
  /** Compute -p. */
  UPolynomial operator-(const UPolynomial& p);
  /** Negate a polynomial in place. */
  void neg(UPolynomial& p);

  /** Add two polynomials. */
  UPolynomial operator+(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Subtract two polynomials. */
  UPolynomial operator-(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Multiply two polynomials. */
  UPolynomial operator*(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Multiply a polynomial and an integer. */
  UPolynomial operator*(const UPolynomial& lhs, const Integer& rhs);
  /** Multiply an integer and a polynomial. */
  UPolynomial operator*(const Integer& lhs, const UPolynomial& rhs);

  /** Compute lhs^rhs. */
  UPolynomial pow(const UPolynomial& lhs, long rhs);

  /** Compute the first derivative. */
  UPolynomial derivative(const UPolynomial& p);

  /** Check whether lhs divides rhs. */
  bool divides(const UPolynomial& lhs, const UPolynomial& rhs);

  /** Divide all degrees within lhs by rhs. All degrees (with non-zero
   * coefficients) must be divisible by rhs. */
  UPolynomial div_degrees(const UPolynomial& lhs, long rhs);

  /** Divides two polynomials, assuming all necessary inverses can be computed.
   */
  UPolynomial div_exact(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Divides a polynomial by an integer, assuming that all coefficients are
   * divisible by rhs. */
  UPolynomial div_exact(const UPolynomial& lhs, const Integer& rhs);
  /** Computes the remainder of two polynomials, assuming all necessary inverses
   * can be computed. */
  UPolynomial rem_exact(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Computes the quotient and the remainder of two polynomials, assuming all
   * necessary inverses can be computed. */
  std::pair<UPolynomial, UPolynomial> div_rem_exact(const UPolynomial& lhs,
                                                    const UPolynomial& rhs);
  /** Performs pseudo-division such that
   *  lc(rhs)^(lhs_deg - rhs_deg + 1) * lhs = div * rhs + rem
   * and returns (div, rem).
   */
  std::pair<UPolynomial, UPolynomial> div_rem_pseudo(const UPolynomial& lhs,
                                                     const UPolynomial& rhs);

  /** Computes the content of a polynomial.
   * The content is the gcd of all coefficients with the sign of the leading
   * coefficient.
   */
  Integer content(const UPolynomial& p);
  /** Makes a polynomial primitive in place.
   * A polynomial is primitive if its content is one.
   */
  void make_primitive(UPolynomial& p);
  /** Computes the primitive part, that is p / content(p). */
  UPolynomial primitive_part(const UPolynomial& p);

  /** Compute the greatest common divisor. */
  UPolynomial gcd(const UPolynomial& lhs, const UPolynomial& rhs);
  /** Compute the extended greatest common divisor. */
  UPolynomial extended_gcd(const UPolynomial& lhs, const UPolynomial& rhs,
                           UPolynomial& u, UPolynomial& v);
  /** Solves the equations
   *  u * p + v * q = r
   * assuming that gcd(p,q) divides r.
   */
  void solve_bezout(const UPolynomial& p, const UPolynomial& q,
                    const UPolynomial& r, UPolynomial& u, UPolynomial& v);

  /** Compute a square-free factorization of a polynomial. */
  std::vector<UPolynomial> square_free_factors(const UPolynomial& p,
                                               bool with_constant = false);

  /** Compute the sturm sequence of a polynomial. */
  std::vector<UPolynomial> sturm_sequence(const UPolynomial& p);

  class AlgebraicNumber;
  class RationalInterval;

  /** Count the real roots of a polynomial within a rational interval.*/
  std::size_t count_real_roots(const UPolynomial& p,
                               const RationalInterval& ri);

  /** Isolate the real roots of a UPolynomial, returning them as algebraic
   * numbers. The roots are sorted in increasing numerical order.
   */
  std::vector<AlgebraicNumber> isolate_real_roots(const UPolynomial& p);

  /**
   * Finds the roots for a polynomial mod p using rabin root finding.
   */
  std::vector<Integer> find_roots_Zp(const UPolynomial& p);

}  // namespace poly
