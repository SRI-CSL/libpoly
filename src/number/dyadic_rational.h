/*
 * dyadic_rational.h
 *
 *  Created on: Jan 22, 2014
 *      Author: dejan
 */

#pragma once

#include <stdio.h>
#include <gmp.h>

#include "integer.h"

/**
 * Fraction of the form a/2^n, with n non-negative. The fraction is reduced, i.e.
 * if n is positive a is not divisible by 2.
 */
typedef struct {
  __mpz_struct a;
  unsigned long n;
} dyadic_rational_t;

/** Interface to the dyadic rationals */
typedef struct {

  /**
   * Construct rational 0/1.
   */
  void (*construct) (dyadic_rational_t* q);

  /**
   * Construct rational a*2^n.
   */
  void (*construct_from_int) (dyadic_rational_t* q, long a, unsigned long n);

  /**
   * Construct rational a/1.
   */
  void (*construct_from_integer) (dyadic_rational_t* q, const integer_t* a);


  /**
   * Construct the dyadic representation of x.
   */
  void (*construct_from_double) (dyadic_rational_t* q, double x);

  /**
   * Construct a copy.
   */
  void (*construct_copy) (dyadic_rational_t* q, const dyadic_rational_t* from);

  /**
   * Assign the given value.
   */
  void (*assign) (dyadic_rational_t* q, const dyadic_rational_t* from);

  /**
   * Assign the rational a given a*2^n
   */
  void (*assign_int) (dyadic_rational_t* q, long a, unsigned long n);

  /**
   * Deallocates the number.
   */
  void (*destruct) (dyadic_rational_t* q);

  /**
   * Prints the number to the given stream.
   */
  int (*print) (const dyadic_rational_t* c, FILE* out);

  /**
   * Returns the string representation of the number.
   */
  char* (*to_string) (const dyadic_rational_t* q);

  /**
   * Return the double representation of the rational.
   */
  double (*to_double) (const dyadic_rational_t* q);

  /**
   * Returns the sign of the rational.
   */
  int (*sgn) (const dyadic_rational_t* q);

  /**
   * Compare the two numbers.
   */
  int (*cmp) (const dyadic_rational_t* q1, const dyadic_rational_t* q2);

  /**
   * Swap two numbers.
   */
  void (*swap) (dyadic_rational_t* q1, dyadic_rational_t* q2);

  /**
   * Compute sum = a + b.
   */
  void (*add) (dyadic_rational_t* sum, const dyadic_rational_t* a, const dyadic_rational_t* b);

  /**
   * Compute sum = a + b.
   */
  void (*add_integer) (dyadic_rational_t* sum, const dyadic_rational_t* a, const integer_t* b);

  /**
   * Compute sub = a - b.
   */
  void (*sub) (dyadic_rational_t* sub, const dyadic_rational_t* a, const dyadic_rational_t* b);

  /**
   * Compute neg = -a.
   */
  void (*neg) (dyadic_rational_t* neg, const dyadic_rational_t* a);

  /**
   * Compute product = a * b.
   */
  void (*mul) (dyadic_rational_t* mul, const dyadic_rational_t* a, const dyadic_rational_t* b);

  /**
   * Compute product = a*2^n
   */
  void (*mul_2exp) (dyadic_rational_t* mul, const dyadic_rational_t* a, unsigned long n);

  /**
   * Compute power = a^n in the given ring.
   */
  void (*pow) (dyadic_rational_t* pow, const dyadic_rational_t*a, unsigned long n);

  /**
   * Compute div = a/2^n
   */
  void (*div_2exp) (dyadic_rational_t* div, const dyadic_rational_t* a, unsigned long n);

} dyadic_rational_ops_t;

/** Implementation */
extern const dyadic_rational_ops_t dyadic_rational_ops;
