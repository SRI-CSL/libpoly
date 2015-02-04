/*
 * rational.h
 *
 *  Created on: Jan 13, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>
#include <gmp.h>

#include "integer.h"
#include "dyadic_rational.h"

/** Use the GMP rationals */
typedef __mpq_struct lp_rational_t;

/** Interface to the rationals */
typedef struct {

  /**
   * Construct rational 0/1.
   */
  void (*construct) (lp_rational_t* q);

  /**
   * Construct rational a/b.
   */
  void (*construct_from_int) (lp_rational_t* q, long a, unsigned long b);

  /**
   * Construct rational a/1.
   */
  void (*construct_from_integer) (lp_rational_t* q, const lp_integer_t* a);

  /**
   * Construct a rational representation of x.
   */
  void (*construct_from_double) (lp_rational_t* q, double x);

  /**
   * Construct a rational number from a dyadic rational.
   */
  void (*construct_from_dyadic) (lp_rational_t* q, const lp_dyadic_rational_t* qd);

  /**
   * Construct a copy of the given coefficient. The coefficient will be
   * normalized according to the given ring.
   */
  void (*construct_copy) (lp_rational_t* q, const lp_rational_t* from);

  /**
   * Assign the rational a given rational.
   */
  void (*assign) (lp_rational_t* q, const lp_rational_t* from);

  /**
   * Assign the rational a given a/b.
   */
  void (*assign_int) (lp_rational_t* q, long a, unsigned long b);

  /**
   * Deallocates the rational.
   */
  void (*destruct) (lp_rational_t* q);

  /**
   * Prints the rational to the given stream.
   */
  int (*print) (const lp_rational_t* c, FILE* out);

  /**
   * Returns the string representation of the rational.
   */
  char* (*to_string) (const lp_rational_t* q);

  /**
   * Return the double representation of the rational.
   */
  double (*to_double) (const lp_rational_t* q);

  /**
   * Returns the sign of the rational.
   */
  int (*sgn) (const lp_rational_t* q);

  /**
   * Compare the two rationals.
   */
  int (*cmp) (const lp_rational_t* q1, const lp_rational_t* q2);

  /**
   * Swap two rationals.
   */
  void (*swap) (lp_rational_t* q1, lp_rational_t* q2);

  /**
   * Compute sum = a + b.
   */
  void (*add) (lp_rational_t* sum, const lp_rational_t* a, const lp_rational_t* b);

  /**
   * Compute sum = a + b.
   */
  void (*add_integer) (lp_rational_t* sum, const lp_rational_t* a, const lp_integer_t* b);

  /**
   * Compute sub = a - b.
   */
  void (*sub) (lp_rational_t* sub, const lp_rational_t* a, const lp_rational_t* b);

  /**
   * Compute neg = -a.
   */
  void (*neg) (lp_rational_t* neg, const lp_rational_t* a);

  /**
   * Compute the inverse of a. Assumes a != 0.
   */
  void (*inv) (lp_rational_t* inv, const lp_rational_t* a);

  /**
   * Compute product = a * b.
   */
  void (*mul) (lp_rational_t* mul, const lp_rational_t* a, const lp_rational_t* b);

  /**
   * Compute product = a*2^n
   */
  void (*mul_2exp) (lp_rational_t* mul, const lp_rational_t* a, unsigned n);

  /**
   * Compute power = a^n in the given ring.
   */
  void (*pow) (lp_rational_t* pow, const lp_rational_t*a, unsigned n);

  /**
   * Compute a/b = div.
   */
  void (*div) (lp_rational_t* div, const lp_rational_t* a, const lp_rational_t* b);

  /**
   * Compute div = a/2^n
   */
  void (*div_2exp) (lp_rational_t* div, const lp_rational_t* a, unsigned n);

} lp_rational_ops_t;

/** Implementation */
extern const lp_rational_ops_t lp_rational_ops;
