/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "poly.h"

#ifdef __cplusplus
extern "C" {
#endif

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

/** Default ring is the whole of integers */
extern lp_int_ring_t* lp_Z;

/**
 * Create a new ring. The new modulus is attached, so in order to remove
 * it you need to detach it. The ring is reference counted, so it will get
 * deallocated when the last user detaches it.
 */
lp_int_ring_t* lp_int_ring_create(const lp_integer_t* M, int is_prime);

/** Attach the ring. */
void lp_int_ring_attach(lp_int_ring_t* K);
/** Detach the ring. */
void lp_int_ring_detach(lp_int_ring_t* K);

/** Are the two rings equal */
int lp_int_ring_equal(const lp_int_ring_t* K1, const lp_int_ring_t* K2);

/** Print */
int lp_int_ring_print(const lp_int_ring_t* K, FILE* out);

/** Get the string representation */
char* lp_int_ring_to_string(const lp_int_ring_t* K);

/** Construct a 0 integer. */
void lp_integer_construct(lp_integer_t* c);

/**
 * Construct a integer from the given rational x. The rational x must be
 * an integer.
 */
void lp_integer_construct_from_rational(const lp_int_ring_t* K, lp_integer_t* c, const lp_rational_t* q);

/**
 * Construct a integer from the given integer x. The integer will be
 * normalized according to the given ring.
 */
void lp_integer_construct_from_int(const lp_int_ring_t* K, lp_integer_t* c, long x);

/**
 * Construct a integer from the given string representation. The
 * integer will be normalized according to the given ring.
 */
void lp_integer_construct_from_string(const lp_int_ring_t* K, lp_integer_t* c, const char* x, int base);

/**
 * Construct a copy of the given integer. The integer will be
 * normalized according to the given ring.
 */
void lp_integer_construct_copy(const lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from);

/**
 * Assign the integer a given integer. The integer will be
 * normalized according to the given ring.
 */
void lp_integer_assign(const lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from);

/**
 * Assign the integer a given integer. The integer will be
 * normalized according to the given ring.
 */
void lp_integer_assign_int(const lp_int_ring_t* K, lp_integer_t* c, long x);

/**
 * Deallocates the integer.
 */
void lp_integer_destruct(lp_integer_t* c);

/**
 * Prints the integer to the given stream.
 */
int lp_integer_print(const lp_integer_t* c, FILE* out);

/**
 * Returns the number of bits needed for this number.
 */
size_t lp_integer_bits(const lp_integer_t* c);

/**
 * Prints the integer matrix (m)x(n) to the given stream.
 */
int lp_integer_print_matrix(const lp_integer_t* c_array, size_t m, size_t n, FILE* out);

/**
 * Returns the string representation of the integer.
 */
char* lp_integer_to_string(const lp_integer_t* c);

/**
 * Returns true if the integer fits in long.
 */
int lp_integer_fits_int(const lp_integer_t* c);

/**
 * Returns the int representation of the integer (truncated but keeps sign).
 */
long lp_integer_to_int(const lp_integer_t* c);

/**
 * Returns the double representation of the integer.
 */
double lp_integer_to_double(const lp_integer_t* c);

/**
 * Returns true if the integer is a prime number.
 */
int lp_integer_is_prime(const lp_integer_t* c);

/**
 * returns true if the integer is zero.
 */
int lp_integer_is_zero(const lp_int_ring_t* K, const lp_integer_t* c);

/**
 * Returns true if the integer is in the given ring by value.
 */
int lp_integer_in_ring(const lp_int_ring_t* K, const lp_integer_t* c);

/**
 * Returns the sign of the integer. The sign is depends on the ring that
 * the integer was created in (the given ring). In a modular ring, the
 * sign is negative if > floor(M/2) when normalized.
 */
int lp_integer_sgn(const lp_int_ring_t* K, const lp_integer_t* c);

/**
 * Compare the two integers in the ring. Not necessarily +/- 1, could be
 * any integer, only the sign matters.
 */
int lp_integer_cmp(const lp_int_ring_t* K, const lp_integer_t* c, const lp_integer_t* to);

/**
 * Compare the two integers in the ring. Same as for cmp.
 */
int lp_integer_cmp_int(const lp_int_ring_t* K, const lp_integer_t* c, long to);

/**
 * Returns true if a divides b in the given ring. In Z this is regular
 * division. In a modular ring it depends on the modulus. In a modular ring
 * with modulus M a divides b iff gcd(a, M) divides b. If the ring is prime
 * any non-zero a divides any b (it's a field).
 */
int lp_integer_divides(const lp_int_ring_t* K, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Swap two integers.
 */
void lp_integer_swap(lp_integer_t* a, lp_integer_t* b);

/**
 * Compute a ++.
 */
void lp_integer_inc(const lp_int_ring_t* K, lp_integer_t* a);

/**
 * Compute a --.
 */
void lp_integer_dec(const lp_int_ring_t* K, lp_integer_t* a);

/**
 * Compute sum = a + b in the given ring.
 */
void lp_integer_add(const lp_int_ring_t* K, lp_integer_t* sum, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute sub = a - b in the given ring.
 */
void lp_integer_sub(const lp_int_ring_t* K, lp_integer_t* sub, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute neg = -a in the given ring. Not that, for example in Z_4, for a =
 * 2, we get that neg = 2. In Z_3 the numbers are in {-1, 0, 1, 2}, and -2 is
 * represented as 2.
 */
void lp_integer_neg(const lp_int_ring_t* K, lp_integer_t* neg, const lp_integer_t* a);

/**
 * Compute the absolute value.
 */
void lp_integer_abs(const lp_int_ring_t* K, lp_integer_t* abs, const lp_integer_t* a);

/**
 * Compute the inverse of a in the given ring. Assumes it has an inverse.
 */
void lp_integer_inv(const lp_int_ring_t* K, lp_integer_t* inv, const lp_integer_t* a);

/**
 * Compute product = a * b in the given ring.
 */
void lp_integer_mul(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute product = a * b in the given ring.
 */
void lp_integer_mul_int(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, long b);

/**
 * Compute product = a*2^n
 */
void lp_integer_mul_pow2(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, unsigned n);

/**
 * Compute power = a^n in the given ring.
 */
void lp_integer_pow(const lp_int_ring_t* K, lp_integer_t* pow, const lp_integer_t* a, unsigned n);

/**
 * Compute the square root of a (in Z).
 */
void lp_integer_sqrt_Z(lp_integer_t* sqrt, const lp_integer_t* a);

/**
 * Compute sum_product += a*b in the given ring.
 */
void lp_integer_add_mul(const lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute sum_product += a*b in the given ring.
 */
void lp_integer_add_mul_int(const lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, int b);

/**
 * Compute sub_product -= a*b in the given ring.
 */
void lp_integer_sub_mul(const lp_int_ring_t* K, lp_integer_t* sub_product, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute a = div*b, in the given ring (assumes that b divides a).
 */
void lp_integer_div_exact(const lp_int_ring_t* K, lp_integer_t* div_Z, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute a = div*b + rem, rounding div towards zero, and r will have the
 * same sign as b, with 0 <= |rem| < |b|.
 */
void lp_integer_div_Z(lp_integer_t* div, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute a = div*b + rem, rounding div towards zero, and r will have the
 * same sign as b, with 0 <= |rem| < |b|.
 */
void lp_integer_rem_Z(lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute a = div*b + rem, rounding div towards zero, and r will have the
 * same sign as b, with 0 <= |rem| < |b|.
 */
void lp_integer_div_rem_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute a = div*2^n + rem, rounding div towards zero, with 0 <= rem < 2^n.
 */
void lp_integer_div_rem_pow2_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, unsigned n);

/**
 * Compute the greatest common divisor (positive) in Z.
 */
void lp_integer_gcd_Z(lp_integer_t* gcd, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Compute the least common multiple (positive) in Z.
 */
void lp_integer_lcm_Z(lp_integer_t* lcm, const lp_integer_t* a, const lp_integer_t* b);

/**
 * Returns the hash of the integer.
 */
size_t lp_integer_hash(const lp_integer_t* a);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
