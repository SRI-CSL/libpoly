/*
 * algebraic_number.h
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>

#include "upolynomial.h"
#include "dyadic_rational.h"

/**
 * Algebraic number represented as the only root of the polynomial f in the
 * interval (a,b). The signs at the points a and b are kept to improve
 * refinement of the interval when needed. If f is 0, then the interval is
 * a single point, and that is the value of the number.
 */
struct lp_algebraic_number_struct {
  lp_upolynomial_t* f;
  lp_dyadic_interval_t I;
  int sgn_at_a, sgn_at_b;
};

/** Operations on algebraic numbers */
typedef struct {

  /**
   * Construct the algebraic number given it's polynomial and the isolating
   * interval. The number takes over the reference of f.
   */
  void (*construct) (lp_algebraic_number_t* a, lp_upolynomial_t* f, const lp_dyadic_interval_t* I);

  /** Construct a zero algebraic number. */
  void (*construct_zero) (lp_algebraic_number_t* a);

  /** Construct a copy of the algebraic number. */
  void (*construct_copy) (lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

  /** Construct the algebraic number from a dyadic rational. */
  void (*construct_from_dyadic_rational) (lp_algebraic_number_t* a, const lp_dyadic_rational_t* q);

  /** Destruct the number. */
  void (*destruct) (lp_algebraic_number_t* a);

  /** Compare two algebraic numbers. */
  int (*cmp) (const lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

  /** Void version of the comparison, use with care. */
  int (*cmp_void) (const void* a1, const void* a2);

  /** Print the number. */
  int (*print) (const lp_algebraic_number_t* a, FILE* out);

  /** Return a string representation of the number. */
  char* (*to_string) (const lp_algebraic_number_t* a);

  /** Convert to double with the given precision. */
  double (*to_double) (const lp_algebraic_number_t* a);

  /** Refine the number by halfing it's interval. */
  void (*refine) (lp_algebraic_number_t* a);

  /** Addition */
  void (*add) (lp_algebraic_number_t* sum, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

  /** Subtraction */
  void (*sub) (lp_algebraic_number_t* sub, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

  /** Negation */
  void (*neg) (lp_algebraic_number_t* neg, const lp_algebraic_number_t* a);

  /** Multiplication */
  void (*mul) (lp_algebraic_number_t* mul, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

  /** Multiplication */
  void (*pow) (lp_algebraic_number_t* pow, const lp_algebraic_number_t* a, unsigned n);

} lp_algebraic_number_ops_t;

extern const lp_algebraic_number_ops_t lp_algebraic_number_ops;
