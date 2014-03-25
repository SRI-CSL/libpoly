/*
 * algebraic_number.h
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#pragma once

#include <stdio.h>

#include "upolynomial.h"
#include "dyadic_rational.h"

/**
 * Algebraic number represented as the only root of the polynomial f in the
 * interval (a,b). The signs at the points a and b are kept to improve
 * refinement of the interval when needed. If f is 0, then the interval is
 * a single point, and that is the value of the number.
 */
struct algebraic_number_struct {
  upolynomial_t* f;
  dyadic_interval_t I;
  int sgn_at_a, sgn_at_b;
};

typedef struct algebraic_number_struct algebraic_number_t;

/** Operations on algebraic numbers */
typedef struct {

  /**
   * Construct the algebraic number given it's polynomial and the isolating
   * interval. The number takes over the reference of f.
   */
  void (*construct) (algebraic_number_t* a, upolynomial_t* f, const dyadic_interval_t* I);

  /**
   * Construct a zero algebraic number.
   */
  void (*construct_zero) (algebraic_number_t* a);


  /**
   * Construct a copy of the algebraic number.
   */
  void (*construct_copy) (algebraic_number_t* a1, const algebraic_number_t* a2);

  /**
   * Construct the algebraic number from a dyadic rational.
   */
  void (*construct_from_dyadic_rational) (algebraic_number_t* a, const dyadic_rational_t* q);

  /**
   * Destruct the number.
   */
  void (*destruct) (algebraic_number_t* a);

  /**
   * Compare two algebraic numbers.
   */
  int (*cmp) (const algebraic_number_t* a1, const algebraic_number_t* a2);

  /**
   * Void version of the comparison, use with care.
   */
  int (*cmp_void) (const void* a1, const void* a2);

  /**
   * Print the number.
   */
  int (*print) (const algebraic_number_t* a, FILE* out);

  /**
   * Return a string representation of the number.
   */
  char* (*to_string) (const algebraic_number_t* a);

  /**
   * Convert to double with the given precision.
   */
  double (*to_double) (const algebraic_number_t* a);

  /**
   * Refine the number by halfing it's interval.
   */
  void (*refine) (algebraic_number_t* a);

  /** Addition */
  void (*add) (algebraic_number_t* sum, const algebraic_number_t* a, const algebraic_number_t* b);

  /** Subtraction */
  void (*sub) (algebraic_number_t* sub, const algebraic_number_t* a, const algebraic_number_t* b);

  /** Multiplication */
  void (*mul) (algebraic_number_t* mul, const algebraic_number_t* a, const algebraic_number_t* b);

  /** Multiplication */
  void (*pow) (algebraic_number_t* pow, const algebraic_number_t* a, unsigned n);

} algebraic_number_ops_t;

extern const algebraic_number_ops_t algebraic_number_ops;
