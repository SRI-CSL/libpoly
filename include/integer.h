/*
 * integer.h
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>

/** The numeric type for the base integers */
typedef __mpz_struct lp_integer_t;

/** Modulus for computing in the congruence ring */
typedef struct {

  /** Reference counter for the ring */
  size_t ref_count;
  /** Is the modulus prime */
  int is_prime;
  /** Modulus */
  lp_integer_t M;
  /** Lower bound -floor(M-1/2)*/
  lp_integer_t lb;
  /** Upper bound floor(M/2)*/
  lp_integer_t ub;

} lp_int_ring_t;

/** The ring type to pass around */
typedef const lp_int_ring_t* lp_int_ring;

/** Default ring is the whole of integers */
extern lp_int_ring lp_Z;

/** Operations on the integer ring */
typedef struct {

  /**
   * Create a new ring. The new modulus is attached, so in order to remove
   * it you need to detach it. The ring is reference counted, so it will get
   * deallocated when the last user detaches it.
   */
  lp_int_ring (*create) (const lp_integer_t* M, int is_prime);

  /** Attach the ring. */
  void (*attach) (lp_int_ring K);
  /** Detach the ring. */
  void (*detach) (lp_int_ring K);

  /** Are the two rings equal */
  int (*equal) (lp_int_ring K1, lp_int_ring K2);

  /** Print */
  int (*print) (lp_int_ring K, FILE* out);
  /** Get the string representation */
  char* (*to_string) (lp_int_ring K);

} lp_int_ring_ops_t;

/** Implementation of the ring operations */
extern lp_int_ring_ops_t lp_int_ring_ops;

/**
 * Interface for integer operations over an arbitrary ring.
 */
typedef struct {

  /**
   * Construct a integer from the given integer x. The integer will be
   * normalized according to the given ring.
   */
  void (*construct_from_int) (lp_int_ring K, lp_integer_t* c, long x);

  /**
   * Construct a integer from the given string representation. The
   * integer will be normalized according to the given ring.
   */
  void (*construct_from_string) (lp_int_ring K, lp_integer_t* c, const char* x, int base);

  /**
   * Construct a copy of the given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*construct_copy) (lp_int_ring K, lp_integer_t* c, const lp_integer_t* from);

  /**
   * Assign the integer a given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*assign) (lp_int_ring K, lp_integer_t* c, const lp_integer_t* from);

  /**
   * Assign the integer a given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*assign_int) (lp_int_ring K, lp_integer_t* c, long x);

  /**
   * Deallocates the integer.
   */
  void (*destruct) (lp_integer_t* c);

  /**
   * Prints the integer to the given stream.
   */
  int (*print) (const lp_integer_t* c, FILE* out);

  /**
   * Returns the number of bits needed for this number.
   */
  size_t (*bits) (const lp_integer_t* c);

  /**
   * Prints the integer matrix (m)x(n) to the given stream.
   */
  int (*print_matrix) (const lp_integer_t* c, size_t m, size_t n, FILE* out);

  /**
   * Returns the string representation of the integer.
   */
  char* (*to_string) (const lp_integer_t* c);

  /**
   * Returns the int representation of the integer.
   */
  long (*to_int) (const lp_integer_t* c);

  /**
   * Returns true if the integer is a prime number.
   */
  int (*is_prime) (const lp_integer_t* c);

  /**
   * returns true if the integer is zero.
   */
  int (*is_zero) (lp_int_ring K, const lp_integer_t* c);

  /**
   * Returns true if the integer is in the given ring by value.
   */
  int (*in_ring) (lp_int_ring K, const lp_integer_t* c);

  /**
   * Returns the sign of the integer. The sign is depends on the ring that
   * the integer was created in (the given ring). In a modular ring, the
   * sign is negative if > floor(M/2) when normalized.
   */
  int (*sgn) (lp_int_ring K, const lp_integer_t* c);

  /**
   * Compare the two integers in the ring. Not necessarily +/- 1, could be
   * any integer, only the sign matters.
   */
  int (*cmp) (lp_int_ring K, const lp_integer_t* c, const lp_integer_t* to);

  /**
   * Compare the two integers in the ring. Same as for sgn.
   */
  int (*cmp_int) (lp_int_ring K, const lp_integer_t* c, long to);

  /**
   * Returns true if a divides b in the given ring. In Z this is regular
   * division. In a modular ring it depends on the modulus. In a modular ring
   * with modulus M a divides b iff gcd(a, M) divides b. If the ring is prime
   * any non-zero a divides any b (it's a field).
   */
  int (*divides) (lp_int_ring K, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Swap two integers.
   */
  void (*swap) (lp_integer_t* a, lp_integer_t* b);

  /**
   * Compute a ++.
   */
  void (*inc) (lp_int_ring K, lp_integer_t* a);

  /**
   * Compute a --.
   */
  void (*dec) (lp_int_ring K, lp_integer_t* a);

  /**
   * Compute sum = a + b in the given ring.
   */
  void (*add) (lp_int_ring K, lp_integer_t* sum, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute sub = a - b in the given ring.
   */
  void (*sub) (lp_int_ring K, lp_integer_t* sub, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute neg = -a in the given ring. Not that, for example in Z_4, for a =
   * 2, we get that neg = 2. In Z_3 the numbers are in {-1, 0, 1, 2}, and -2 is
   * represented as 2.
   */
  void (*neg) (lp_int_ring K, lp_integer_t* neg, const lp_integer_t* a);

  /**
   * Compute the absolute value.
   */
  void (*abs) (lp_int_ring K, lp_integer_t* abs, const lp_integer_t* a);

  /**
   * Compute the inverse of a in the given ring. Assumes it has an inverse.
   */
  void (*inv) (lp_int_ring K, lp_integer_t* inv, const lp_integer_t* a);

  /**
   * Compute product = a * b in the given ring.
   */
  void (*mul) (lp_int_ring K, lp_integer_t* product, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute product = a * b in the given ring.
   */
  void (*mul_int) (lp_int_ring K, lp_integer_t* product, const lp_integer_t* a, long b);

  /**
   * Compute product = a*2^n
   */
  void (*mul_pow2) (lp_int_ring K, lp_integer_t* product, const lp_integer_t* a, unsigned n);

  /**
   * Compute power = a^n in the given ring.
   */
  void (*pow) (lp_int_ring K, lp_integer_t* pow, const lp_integer_t*a, unsigned n);

  /**
   * Compute the square root of a (in Z).
   */
  void (*sqrt_Z) (lp_integer_t* sqrt, const lp_integer_t* a);

  /**
   * Compute sum_product += a*b in the given ring.
   */
  void (*add_mul) (lp_int_ring K, lp_integer_t* sum_product, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute sum_product += a*b in the given ring.
   */
  void (*add_mul_int) (lp_int_ring K, lp_integer_t* sum_product, const lp_integer_t* a, int b);

  /**
   * Compute sub_product -= a*b in the given ring.
   */
  void (*sub_mul) (lp_int_ring K, lp_integer_t* sub_product, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute a = div*b, in the given ring (assumes that b divides a).
   */
  void (*div_exact) (lp_int_ring K, lp_integer_t* div_Z, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*div_Z) (lp_integer_t* div, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*rem_Z) (lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*div_rem_Z) (lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute a = div*2^n + rem, rounding div towards zero, with 0 <= rem < 2^n.
   */
  void (*div_rem_pow2_Z) (lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, unsigned n);

  /**
   * Compute the greatest common divisor (positive) in Z.
   */
  void (*gcd_Z) (lp_integer_t* gcd, const lp_integer_t* a, const lp_integer_t* b);

  /**
   * Compute the least common multiple (positive) in Z.
   */
  void (*lcm_Z) (lp_integer_t* lcm, const lp_integer_t* a, const lp_integer_t* b);

} lp_integer_ops_struct;

/** Implementation: function table for integer operations */
extern const lp_integer_ops_struct lp_integer_ops;
