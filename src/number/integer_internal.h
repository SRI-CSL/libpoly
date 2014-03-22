/*
 * integer_internal.h
 *
 *  Created on: Mar 1, 2014
 *      Author: dejan
 */

#pragma once

#include "integer.h"

#include <assert.h>
#include <malloc.h>
#include <printf.h>

#define __unused(x) ((void)x)

static inline
int integer_in_ring(int_ring K, const integer_t* c) {
  if (K) {
    int sgn = mpz_sgn(c);
    if (sgn == 0) return 1;
    if (sgn > 0 && mpz_cmp(c, &K->ub) <= 0) return 1;
    if (sgn < 0 && mpz_cmp(&K->lb, c) <= 0) return 1;
    return 0;
  } else {
    // Everything is in Z
    return 1;
  }
}

inline static
void integer_ring_normalize(int_ring K, integer_t* c) {
  if (K && !integer_in_ring(K, c)) {
    // Remainder
    integer_t tmp;
    mpz_init(&tmp);
    // c = M*div + tmp, with 0 < |tmp| < M tmp same sign as c
    mpz_tdiv_r(&tmp, c, &K->M);
    // Swap rem and c
    mpz_swap(c, &tmp);
    // Get the sign of c
    int sgn = mpz_sgn(c);
    // Make smaller than the upper bound
    if (sgn > 0 && mpz_cmp(c, &K->ub) > 0) {
      mpz_sub(&tmp, c, &K->M);
      mpz_swap(c, &tmp);
    }
    // Make bigger than the lower bound
    if (sgn < 0 && mpz_cmp(c, &K->lb) < 0) {
      // For negative ones, we might have to subtract
      mpz_add(&tmp, c, &K->M);
      mpz_swap(c, &tmp);
    }
    // Remove the temp
    mpz_clear(&tmp);
    assert(integer_in_ring(K, c));
  }
}

static inline
void integer_construct_from_int(int_ring K, integer_t* c, long x) {
  mpz_init_set_si(c, x);
  integer_ring_normalize(K, c);
}

static inline
void integer_construct_copy(int_ring K, integer_t* c, const integer_t* from) {
  mpz_init_set(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_assign(int_ring K, integer_t* c, const integer_t* from) {
  mpz_set(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_assign_int(int_ring K, integer_t* c, long from) {
  mpz_set_si(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_destruct(integer_t* c) {
  mpz_clear(c);
}

static inline
int integer_print(const integer_t* c, FILE* out) {
  return mpz_out_str(out, 10, c);
}

static inline
int integer_bits(const integer_t* c) {
  return mpz_sizeinbase(c, 2);
}

static inline
long integer_to_int(const integer_t* c) {
  return mpz_get_si(c);
}

static inline
int integer_is_prime(const integer_t* c) {
  return mpz_probab_prime_p(c, 25);
}

static inline
int integer_is_zero(int_ring K, const integer_t* c) {
  if (K) {
    integer_t c_normalized;
    integer_construct_copy(K, &c_normalized, c);
    int sgn = mpz_sgn(&c_normalized);
    integer_destruct(&c_normalized);
    return sgn == 0;
  } else {
    return mpz_sgn(c) == 0;
  }
}

static inline
int integer_sgn(int_ring K, const integer_t* c) {
  if (K) {
    integer_t c_normalized;
    integer_construct_copy(K, &c_normalized, c);
    int sgn = mpz_sgn(&c_normalized);
    integer_destruct(&c_normalized);
    return sgn;
  } else {
    return mpz_sgn(c);
  }
}

static inline
int integer_cmp(int_ring K, const integer_t* c, const integer_t* to) {
  if (K) {
    integer_t c_normalized, to_normalized;
    integer_construct_copy(K, &c_normalized, c);
    integer_construct_copy(K, &to_normalized, to);
    int cmp = mpz_cmp(&c_normalized, &to_normalized);
    integer_destruct(&c_normalized);
    integer_destruct(&to_normalized);
    return cmp;
  } else {
    return mpz_cmp(c, to);
  }
}

static inline
int integer_cmp_int(int_ring K, const integer_t* c, long to) {
  if (K) {
    integer_t c_normalized, to_normalized;
    integer_construct_copy(K, &c_normalized, c);
    integer_construct_from_int(K, &to_normalized, to);
    int cmp = mpz_cmp(&c_normalized, &to_normalized);
    integer_destruct(&c_normalized);
    integer_destruct(&to_normalized);
    return cmp;
  } else {
    return mpz_cmp_si(c, to);
  }
}

static inline
int integer_divides(int_ring K, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  if (K) {
    // In a prime ring, it's always divisible
    if (K->is_prime) return integer_ops.sgn(Z, a);
    // Otherwise compute the gcd
    integer_t gcd;
    mpz_init(&gcd);
    mpz_gcd(&gcd, a, &K->M);
    int divides = mpz_divisible_p(b, &gcd);
    mpz_clear(&gcd);
    return divides;
  } else {
    return mpz_divisible_p(b, a);
  }
}

static inline
void integer_swap(int_ring K, integer_t* a, integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_swap(a, b);
}

static inline
void integer_inc(int_ring K, integer_t* a) {
  assert(integer_in_ring(K, a));
  integer_t tmp;
  mpz_init(&tmp);
  mpz_add_ui(&tmp, a, 1);
  mpz_swap(&tmp, a);
  mpz_clear(&tmp);
  integer_ring_normalize(K, a);
}

static inline
void integer_dec(int_ring K, integer_t* a) {
  assert(integer_in_ring(K, a));
  integer_t tmp;
  mpz_init(&tmp);
  mpz_sub_ui(&tmp, a, 1);
  mpz_swap(&tmp, a);
  mpz_clear(&tmp);
  integer_ring_normalize(K, a);
}

static inline
void integer_add(int_ring K, integer_t* sum, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_add(sum, a, b);
  integer_ring_normalize(K, sum);
}

static inline
void integer_sub(int_ring K, integer_t* sub, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_sub(sub, a, b);
  integer_ring_normalize(K, sub);
}

static inline
void integer_neg(int_ring K, integer_t* neg, const integer_t* a) {
  assert(integer_in_ring(K, a));
  mpz_neg(neg, a);
  integer_ring_normalize(K, neg);
}

static inline
void integer_abs(int_ring K, integer_t* abs, const integer_t* a) {
  assert(integer_in_ring(K, a));
  mpz_abs(abs, a);
  integer_ring_normalize(K, abs);
}

static inline
void integer_inv(int_ring K, integer_t* inv, const integer_t* a) {
  assert(K);
  assert(integer_in_ring(K, a));
  int result = mpz_invert(inv, a, &K->M);
  assert(result);
  __unused(result);
  integer_ring_normalize(K, inv);
}

static inline
void integer_mul(int_ring K, integer_t* product, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_mul(product, a, b);
  integer_ring_normalize(K, product);
}

static inline
void integer_mul_int(int_ring K, integer_t* product, const integer_t* a, long b) {
  assert(integer_in_ring(K, a));
  mpz_mul_si(product, a, b);
  integer_ring_normalize(K, product);
}

static inline
void integer_mul_pow2(int_ring K, integer_t* power, const integer_t* a, unsigned n) {
  assert(integer_in_ring(K, a));
  assert(n > 0);
  mpz_mul_2exp(power, a, n);
  integer_ring_normalize(K, power);
}

static inline
void integer_pow(int_ring K, integer_t* power, const integer_t*a, unsigned n) {
  assert(integer_in_ring(K, a));
  if (K) {
    mpz_powm_ui(power, a, n, &K->M);
    integer_ring_normalize(K, power);
  } else {
    mpz_pow_ui(power, a, n);
  }
}

static inline
void integer_sqrt_Z(integer_t* sqrt, const integer_t* a) {
  mpz_sqrt(sqrt, a);
}

static inline
void integer_add_mul(int_ring K, integer_t* sum_product, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, sum_product) && integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_addmul(sum_product, a, b);
  integer_ring_normalize(K, sum_product);
}

static inline
void integer_sub_mul(int_ring K, integer_t* sub_product, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, sub_product) && integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_submul(sub_product, a, b);
  integer_ring_normalize(K, sub_product);
}

static inline
void integer_add_mul_int(int_ring K, integer_t* sum_product, const integer_t* a, int b) {
  assert(integer_in_ring(K, sum_product));
  assert(integer_in_ring(K, a));
  if (b > 0) {
    mpz_addmul_ui(sum_product, a, b);
  } else {
    mpz_submul_ui(sum_product, a, -b);
  }
  integer_ring_normalize(K, sum_product);
}

static inline
void integer_div_exact(int_ring K, integer_t* div, const integer_t* a, const integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  if (K) {
    // Solving a = div*b (mod M). Let d = gcd(b, M) with extended gcd, we have
    // that c1*b+c2*M = d. Since d should divide a, the we get the solution
    // multiplying by a/d, obtaining (c1*a/d)*b = a.
    integer_t c1, c2, gcd;
    mpz_init(&c1); mpz_init(&c2); mpz_init(&gcd);
    mpz_gcdext(&gcd, &c1, &c2, b, &K->M);
    assert(mpz_divisible_p(a, &gcd));
    mpz_divexact(&c2, a, &gcd);
    mpz_mul(div, &c1, &c2);
    mpz_clear(&c1); mpz_clear(&c2); mpz_clear(&gcd);
    integer_ring_normalize(K, div);
  } else {
    mpz_divexact(div, a, b);
  }
}

static inline
void integer_div_Z(integer_t* div, const integer_t* a, const integer_t* b) {
  mpz_tdiv_q(div, a, b);
}

static inline
void integer_rem_Z(integer_t* rem, const integer_t* a, const integer_t* b) {
  mpz_tdiv_r(rem, a, b);
}

static inline
void integer_div_rem_Z(integer_t* div, integer_t* rem, const integer_t* a, const integer_t* b) {
  mpz_tdiv_qr(div, rem, a, b);
}

static inline
void integer_div_rem_pow2_Z(integer_t* div, integer_t* rem, const integer_t* a, unsigned n) {
  mpz_tdiv_q_2exp(div, a, n);
  mpz_tdiv_r_2exp(rem, a, n);
}

static inline
void integer_gcd_Z(integer_t* gcd, const integer_t* a, const integer_t* b) {
  mpz_gcd(gcd, a, b);
}

static inline
void integer_lcm_Z(integer_t* lcm, const integer_t* a, const integer_t* b) {
  mpz_lcm(lcm, a, b);
}
