#pragma once

#include <iosfwd>
#include <vector>

#include "../polynomial.h"
#include "assignment.h"
#include "integer.h"
#include "interval.h"
#include "interval_assignment.h"
#include "sign_condition.h"
#include "value.h"
#include "variable.h"

namespace poly {

  /**
   * Implements a wrapper for lp_polynomial_t.
   */
  class Polynomial {
    /** The actual polynomial. */
    deleting_unique_ptr<lp_polynomial_t> mPoly;

   public:
    /** Create from a lp_polynomial_t pointer, claiming it's ownership. */
    Polynomial(lp_polynomial_t* poly);
    /** Copy from an internal lp_polynomial_t pointer. */
    Polynomial(const lp_polynomial_t* poly);
    /** Construct a zero polynomial from an internal context pointer. */
    Polynomial(const lp_polynomial_context_t* c);
    /** Construct a zero polynomial from a custom context. */
    Polynomial(const Context& c);
    /** Construct a zero polynomial. */
    Polynomial();

    /** Construct from a variable and a custom context. */
    Polynomial(const Context& c, Variable v);
    /** Construct from a variable. */
    Polynomial(Variable v);

    /** Construct i * v^n from a custom context. */
    Polynomial(const Context& c, const Integer& i, Variable v, unsigned n);
    /** Construct i * v^n. */
    Polynomial(const Integer& i, Variable v, unsigned n);

    /** Construct from an integer and a custom context. */
    Polynomial(const Context& c, const Integer& i);
    /** Construct from an integer. */
    Polynomial(const Integer& i);

    /** Construct from an integer and a custom context. */
    Polynomial(const Context& c, long i);
    /** Construct from an integer. */
    Polynomial(long i);

    /** Copy from a Polynomial. */
    Polynomial(const Polynomial& p);
    /** Move from a Polynomial. */
    Polynomial(Polynomial&& p) noexcept;

    /** Copy from a Polynomial. */
    Polynomial& operator=(const Polynomial& p);
    /** Move from a Polynomial. */
    Polynomial& operator=(Polynomial&& p) noexcept;

    /** Get a non-const pointer to the internal lp_polynomial_t. Handle with
     * care!
     */
    lp_polynomial_t* get_internal();
    /** Get a const pointer to the internal lp_polynomial_t. */
    const lp_polynomial_t* get_internal() const;
    /** Release the lp_polynomial_t pointer. This yields ownership of the
     * returned pointer. */
    lp_polynomial_t* release();
  };

  /** Swap two polynomials. */
  void swap(Polynomial& lhs, Polynomial& rhs);

  /** Return the hash of a polynomial. */
  std::size_t hash(const Polynomial& p);

  /** Stream the given Polynomial to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Polynomial& p);

  /** Check if the given polynomial is zero. */
  bool is_zero(const Polynomial& p);
  /** Check if the given polynomial is constant. */
  bool is_constant(const Polynomial& p);
  /** Check if the given polynomial is linear. */
  bool is_linear(const Polynomial& p);
  /** Check if the leading coefficient is constant. */
  bool is_lc_constant(const Polynomial& p);
  /** Get the sign of the leading coefficient. Assumes is_lc_constant(p). */
  int lc_sgn(const Polynomial& p);
  /** Obtain the degree of the given polynomial in its main variable. */
  std::size_t degree(const Polynomial& p);
  /** Obtain the main variable of the given polynomial. */
  Variable main_variable(const Polynomial& p);
  /** Obtain the k'th coefficient of a polynomial. */
  Polynomial coefficient(const Polynomial& p, std::size_t k);
  /** Obtain the leading coefficient of a polynomial. */
  Polynomial leading_coefficient(const Polynomial& p);
  /** Obtain all non-constant coefficients of a polynomial. */
  std::vector<Polynomial> coefficients(const Polynomial& p);

  /** Check if the given polynomial is univariate. */
  bool is_univariate(const Polynomial& p);
  /** Converts a polynomial to a univariate polynomial. Assumes
   * is_univariate(p). */
  UPolynomial to_univariate(const Polynomial& p);
  /** Check if the given polynomial is univariate over a given assignment and
   * the main variable is unassigned. */
  bool is_univariate_over_assignment(const Polynomial& p, const Assignment& a);
  /** Check if the given polynomial is completely assigned over a given
   * assignment. */
  bool is_assigned_over_assignment(const Polynomial& p, const Assignment& a);
  /** Evaluates p over a given assignment and returns an univariate polynomial.
   * Assumes that a assigns all variable in p but the top variable.
   * Assumes that a assigns to integer only. */
  UPolynomial to_univariate(const Polynomial& p, const Assignment& a);
  /** Compute the sign of a polynomial over an assignment. */
  int sgn(const Polynomial& p, const Assignment& a);
  /** Evaluate a polynomial over an assignment. */
  Value evaluate(const Polynomial& p, const Assignment& a);
  /** Evaluate a polynomial constraint over an assignment. */
  bool evaluate_constraint(const Polynomial& p, const Assignment& a,
                           SignCondition sc);
  /** Evaluate a polynomial over an interval assignment. The result is only an
   * approximation. */
  Interval evaluate(const Polynomial& p, const IntervalAssignment& a);

  /** Compare polynomials. */
  bool operator==(const Polynomial& lhs, const Polynomial& rhs);
  /** Compare polynomials. */
  bool operator!=(const Polynomial& lhs, const Polynomial& rhs);
  /** Compare polynomials. */
  bool operator<(const Polynomial& lhs, const Polynomial& rhs);
  /** Compare polynomials. */
  bool operator<=(const Polynomial& lhs, const Polynomial& rhs);
  /** Compare polynomials. */
  bool operator>(const Polynomial& lhs, const Polynomial& rhs);
  /** Compare polynomials. */
  bool operator>=(const Polynomial& lhs, const Polynomial& rhs);

  /** Add two polynomials. */
  Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs);
  /** Add a polynomial and an integer. */
  Polynomial operator+(const Polynomial& lhs, const Integer& rhs);
  /** Add an integer and a polynomial. */
  Polynomial operator+(const Integer& lhs, const Polynomial& rhs);
  /** Add and assign two polynomials. */
  Polynomial& operator+=(Polynomial& lhs, const Polynomial& rhs);
  /** Compute lhs += rhs1 * rhs2. */
  Polynomial& add_mul(Polynomial& lhs, const Polynomial& rhs1, const Polynomial& rhs2);

  /** Unary negation of a polynomial. */
  Polynomial operator-(const Polynomial& p);
  /** Subtract two polynomials. */
  Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs);
  /** Subtract an integer from a polynomial. */
  Polynomial operator-(const Polynomial& lhs, const Integer& rhs);
  /** Subtract a polynomial from an integer. */
  Polynomial operator-(const Integer& lhs, const Polynomial& rhs);
  /** Subtract and assign two polynomials. */
  Polynomial& operator-=(Polynomial& lhs, const Polynomial& rhs);
  /** Compute lhs -= rhs1 * rhs2. */
  Polynomial& sub_mul(Polynomial& lhs, const Polynomial& rhs1, const Polynomial& rhs2);

  /** Multiply two polynomials. */
  Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs);
  /** Multiply a polynomial and an integer. */
  Polynomial operator*(const Polynomial& lhs, const Integer& rhs);
  /** Multiply an integer and a polynomial. */
  Polynomial operator*(const Integer& lhs, const Polynomial& rhs);
  /** Multiply and assign two polynomials. */
  Polynomial& operator*=(Polynomial& lhs, const Polynomial& rhs);

  /** Multiply with x^n where x is the main variable. */
  Polynomial shl(const Polynomial& lhs, unsigned n);
  /** Compute a polynomial to some power. */
  Polynomial pow(const Polynomial& lhs, unsigned exp);

  /** Checks whether lhs divides rhs. */
  bool divides(const Polynomial& lhs, const Polynomial& rhs);
  /** Compute the quotient of two polynomials, assuming divides(rhs, lhs). */
  Polynomial div(const Polynomial& lhs, const Polynomial& rhs);
  /** Compute the remainder of two polynomials. */
  Polynomial rem(const Polynomial& lhs, const Polynomial& rhs);
  /** Compute the pseudo-remainder of two polynomials. */
  Polynomial prem(const Polynomial& lhs, const Polynomial& rhs);
  /** Compute the sparse pseudo-remainder of two polynomials. */
  Polynomial sprem(const Polynomial& lhs, const Polynomial& rhs);
  /** Compute quotient and remainder of two polynomials. */
  std::pair<Polynomial, Polynomial> div_rem(const Polynomial& lhs,
                                            const Polynomial& rhs);

  /** Compute the derivative of a polynomial (in its main variable). */
  Polynomial derivative(const Polynomial& p);

  /** Compute the gcd of two polynomials. */
  Polynomial gcd(const Polynomial& p, const Polynomial& q);
  /** Compute the lcm of two polynomials. */
  Polynomial lcm(const Polynomial& p, const Polynomial& q);

  /** Compute the content of a polynomial. */
  Polynomial content(const Polynomial& p);
  /** Compute the primitive part of a polynomial. */
  Polynomial primitive_part(const Polynomial& p);
  /** Compute the content and the primitive part of a polynomial. */
  std::pair<Polynomial, Polynomial> content_primitive_part(const Polynomial& p);

  /** Compute the resultant of two polynomials. */
  Polynomial resultant(const Polynomial& p, const Polynomial& q);
  /** Compute the discriminant of a polynomial. */
  Polynomial discriminant(const Polynomial& p);

  /** Compute the principal subresultant coefficients of two polynomials. */
  std::vector<Polynomial> psc(const Polynomial& p, const Polynomial& q);

  /**
   * Compute a square-free factorization of a polynomial.
   * Attention: this does not yield a full factorization!
   */
  std::vector<Polynomial> square_free_factors(const Polynomial& p);
  /** Compute a content-free factorization of a polynomial. */
  std::vector<Polynomial> content_free_factors(const Polynomial& p);

  /** Isolate the real roots of a Polynomial with respect to an Assignment for
   * all but the main variable. */
  std::vector<Value> isolate_real_roots(const Polynomial& p,
                                        const Assignment& a);

  std::vector<Interval> infeasible_regions(const Polynomial& p,
                                           const Assignment& a,
                                           SignCondition sc);

}  // namespace poly
