/*
 * rational.h
 *
 *  Created on: Jan 13, 2014
 *      Author: dejan
 */

#pragma once

#include <stdio.h>
#include <gmp.h>

#include "integer.h"
#include "dyadic_rational.h"

/** Use the GMP rationals */
typedef __mpq_struct rational_t;

/** Interface to the rationals */
typedef struct {

  /**
   * Construct rational 0/1.
   */
  void (*construct) (rational_t* q);

  /**
   * Construct rational a/b.
   */
  void (*construct_from_int) (rational_t* q, long a, unsigned long b);

  /**
   * Construct rational a/1.
   */
  void (*construct_from_integer) (rational_t* q, const integer_t* a);

  /**
   * Construct a rational representation of x.
   */
  void (*construct_from_double) (rational_t* q, double x);

  /**
   * Construct a rational number from a dyadic rational.
   */
  void (*construct_from_dyadic) (rational_t* q, const dyadic_rational_t* qd);

  /**
   * Construct a copy of the given coefficient. The coefficient will be
   * normalized according to the given ring.
   */
  void (*construct_copy) (rational_t* q, const rational_t* from);

  /**
   * Assign the rational a given rational.
   */
  void (*assign) (rational_t* q, const rational_t* from);

  /**
   * Assign the rational a given a/b.
   */
  void (*assign_int) (rational_t* q, long a, unsigned long b);

  /**
   * Deallocates the rational.
   */
  void (*destruct) (rational_t* q);

  /**
   * Prints the rational to the given stream.
   */
  int (*print) (const rational_t* c, FILE* out);

  /**
   * Returns the string representation of the rational.
   */
  char* (*to_string) (const rational_t* q);

  /**
   * Return the double representation of the rational.
   */
  double (*to_double) (const rational_t* q);

  /**
   * Returns the sign of the rational.
   */
  int (*sgn) (const rational_t* q);

  /**
   * Compare the two rationals.
   */
  int (*cmp) (const rational_t* q1, const rational_t* q2);

  /**
   * Swap two rationals.
   */
  void (*swap) (rational_t* q1, rational_t* q2);

  /**
   * Compute sum = a + b.
   */
  void (*add) (rational_t* sum, const rational_t* a, const rational_t* b);

  /**
   * Compute sum = a + b.
   */
  void (*add_integer) (rational_t* sum, const rational_t* a, const integer_t* b);

  /**
   * Compute sub = a - b.
   */
  void (*sub) (rational_t* sub, const rational_t* a, const rational_t* b);

  /**
   * Compute neg = -a.
   */
  void (*neg) (rational_t* neg, const rational_t* a);

  /**
   * Compute the inverse of a. Assumes a != 0.
   */
  void (*inv) (rational_t* inv, const rational_t* a);

  /**
   * Compute product = a * b.
   */
  void (*mul) (rational_t* mul, const rational_t* a, const rational_t* b);

  /**
   * Compute product = a*2^n
   */
  void (*mul_2exp) (rational_t* mul, const rational_t* a, unsigned n);

  /**
   * Compute power = a^n in the given ring.
   */
  void (*pow) (rational_t* pow, const rational_t*a, unsigned n);

  /**
   * Compute a/b = div.
   */
  void (*div) (rational_t* div, const rational_t* a, const rational_t* b);

  /**
   * Compute div = a/2^n
   */
  void (*div_2exp) (rational_t* div, const rational_t* a, unsigned n);

} rational_ops_t;

/** Implementation */
extern const rational_ops_t rational_ops;
