/*
 * integer.h
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>

/** The numeric type for the base integers */
typedef __mpz_struct integer_t;

/** Modulus for computing in the congruence ring */
typedef struct {

  /** Reference counter for the ring */
  size_t ref_count;
  /** Is the modulus prime */
  int is_prime;
  /** Modulus */
  integer_t M;
  /** Lower bound -floor(M-1/2)*/
  integer_t lb;
  /** Upper bound floor(M/2)*/
  integer_t ub;

} int_ring_t;

/** The ring type to pass around */
typedef const int_ring_t* int_ring;

/** Default ring is the whole of integers */
int_ring Z;

/** Operations on the integer ring */
typedef struct {

  /**
   * Create a new ring. The new modulus is attached, so in order to remove
   * it you need to detach it. The ring is reference counted, so it will get
   * deallocated when the last user detaches it.
   */
  int_ring (*create) (const integer_t* M, int is_prime);

  /** Attach the ring. */
  void (*attach) (int_ring K);
  /** Detach the ring. */
  void (*detach) (int_ring K);

  /** Are the two rings equal */
  int (*equal) (int_ring K1, int_ring K2);

  /** Print */
  int (*print) (int_ring K, FILE* out);
  /** Get the string representation */
  char* (*to_string) (int_ring K);

} int_ring_ops_t;

/** Implementation of the ring operations */
extern int_ring_ops_t int_ring_ops;

/**
 * Interface for integer operations over an arbitrary ring.
 */
typedef struct {

  /**
   * Construct a integer from the given integer x. The integer will be
   * normalized according to the given ring.
   */
  void (*construct_from_int) (int_ring K, integer_t* c, long x);

  /**
   * Construct a integer from the given string representation. The
   * integer will be normalized according to the given ring.
   */
  void (*construct_from_string) (int_ring K, integer_t* c, const char* x, int base);

  /**
   * Construct a copy of the given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*construct_copy) (int_ring K, integer_t* c, const integer_t* from);

  /**
   * Assign the integer a given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*assign) (int_ring K, integer_t* c, const integer_t* from);

  /**
   * Assign the integer a given integer. The integer will be
   * normalized according to the given ring.
   */
  void (*assign_int) (int_ring K, integer_t* c, long x);

  /**
   * Deallocates the integer.
   */
  void (*destruct) (integer_t* c);

  /**
   * Prints the integer to the given stream.
   */
  int (*print) (const integer_t* c, FILE* out);

  /**
   * Returns the number of bits needed for this number.
   */
  int (*bits) (const integer_t* c);

  /**
   * Prints the integer matrix (m)x(n) to the given stream.
   */
  int (*print_matrix) (const integer_t* c, size_t m, size_t n, FILE* out);

  /**
   * Returns the string representation of the integer.
   */
  char* (*to_string) (const integer_t* c);

  /**
   * Returns the int representation of the integer.
   */
  long (*to_int) (const integer_t* c);

  /**
   * Returns true if the integer is a prime number.
   */
  int (*is_prime) (const integer_t* c);

  /**
   * returns true if the integer is zero.
   */
  int (*is_zero) (int_ring K, const integer_t* c);

  /**
   * Returns true if the integer is in the given ring by value.
   */
  int (*in_ring) (int_ring K, const integer_t* c);

  /**
   * Returns the sign of the integer. The sign is depends on the ring that
   * the integer was created in (the given ring). In a modular ring, the
   * sign is negative if > floor(M/2) when normalized.
   */
  int (*sgn) (int_ring K, const integer_t* c);

  /**
   * Compare the two integers in the ring. Not necessarily +/- 1, could be
   * any integer, only the sign matters.
   */
  int (*cmp) (int_ring K, const integer_t* c, const integer_t* to);

  /**
   * Compare the two integers in the ring. Same as for sgn.
   */
  int (*cmp_int) (int_ring K, const integer_t* c, long to);

  /**
   * Returns true if a divides b in the given ring. In Z this is regular
   * division. In a modular ring it depends on the modulus. In a modular ring
   * with modulus M a divides b iff gcd(a, M) divides b. If the ring is prime
   * any non-zero a divides any b (it's a field).
   */
  int (*divides) (int_ring K, const integer_t* a, const integer_t* b);

  /**
   * Swap two integers. The integers are assumed to already be in the
   * given ring.
   */
  void (*swap) (int_ring K, integer_t* a, integer_t* b);

  /**
   * Compute a ++.
   */
  void (*inc) (int_ring K, integer_t* a);

  /**
   * Compute a --.
   */
  void (*dec) (int_ring K, integer_t* a);

  /**
   * Compute sum = a + b in the given ring.
   */
  void (*add) (int_ring K, integer_t* sum, const integer_t* a, const integer_t* b);

  /**
   * Compute sub = a - b in the given ring.
   */
  void (*sub) (int_ring K, integer_t* sub, const integer_t* a, const integer_t* b);

  /**
   * Compute neg = -a in the given ring. Not that, for example in Z_4, for a =
   * 2, we get that neg = 2. In Z_3 the numbers are in {-1, 0, 1, 2}, and -2 is
   * represented as 2.
   */
  void (*neg) (int_ring K, integer_t* neg, const integer_t* a);

  /**
   * Compute the absolute value.
   */
  void (*abs) (int_ring K, integer_t* abs, const integer_t* a);

  /**
   * Compute the inverse of a in the given ring. Assumes it has an inverse.
   */
  void (*inv) (int_ring K, integer_t* inv, const integer_t* a);

  /**
   * Compute product = a * b in the given ring.
   */
  void (*mul) (int_ring K, integer_t* product, const integer_t* a, const integer_t* b);

  /**
   * Compute product = a * b in the given ring.
   */
  void (*mul_int) (int_ring K, integer_t* product, const integer_t* a, long b);

  /**
   * Compute product = a*2^n
   */
  void (*mul_pow2) (int_ring K, integer_t* product, const integer_t* a, unsigned n);

  /**
   * Compute power = a^n in the given ring.
   */
  void (*pow) (int_ring K, integer_t* pow, const integer_t*a, unsigned n);

  /**
   * Compute the square root of a (in Z).
   */
  void (*sqrt_Z) (integer_t* sqrt, const integer_t* a);

  /**
   * Compute sum_product += a*b in the given ring.
   */
  void (*add_mul) (int_ring K, integer_t* sum_product, const integer_t* a, const integer_t* b);

  /**
   * Compute sum_product += a*b in the given ring.
   */
  void (*add_mul_int) (int_ring K, integer_t* sum_product, const integer_t* a, int b);

  /**
   * Compute sub_product -= a*b in the given ring.
   */
  void (*sub_mul) (int_ring K, integer_t* sub_product, const integer_t* a, const integer_t* b);

  /**
   * Compute a = div*b, in the given ring (assumes that b divides a).
   */
  void (*div_exact) (int_ring K, integer_t* div_Z, const integer_t* a, const integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*div_Z) (integer_t* div, const integer_t* a, const integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*rem_Z) (integer_t* rem, const integer_t* a, const integer_t* b);

  /**
   * Compute a = div*b + rem, rounding div towards zero, and r will have the
   * same sign as b, with 0 <= |rem| < |b|.
   */
  void (*div_rem_Z) (integer_t* div, integer_t* rem, const integer_t* a, const integer_t* b);

  /**
   * Compute a = div*2^n + rem, rounding div towards zero, with 0 <= rem < 2^n.
   */
  void (*div_rem_pow2_Z) (integer_t* div, integer_t* rem, const integer_t* a, unsigned n);

  /**
   * Compute the greatest common divisor (positive) in Z.
   */
  void (*gcd_Z) (integer_t* gcd, const integer_t* a, const integer_t* b);

  /**
   * Compute the least common multiple (positive) in Z.
   */
  void (*lcm_Z) (integer_t* lcm, const integer_t* a, const integer_t* b);

  /**
   * Register printf extension %C for integers.
   */
  void (*register_printf_extension) (void);

} integer_ops_struct;

/** Implementation: function table for integer operations */
extern const integer_ops_struct integer_ops;
