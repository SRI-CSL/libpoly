/*
 * integer_internal.h
 *
 *  Created on: Mar 1, 2014
 *      Author: dejan
 */

#pragma once

#include <integer.h>

#include <assert.h>
#include <malloc.h>
#include <printf.h>

#define __unused(x) ((void)x)

static inline
int integer_in_ring(lp_int_ring_t* K, const lp_integer_t* c) {
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
void integer_ring_normalize(lp_int_ring_t* K, lp_integer_t* c) {
  if (K && !integer_in_ring(K, c)) {
    // Remainder
    lp_integer_t tmp;
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
void integer_construct_from_int(lp_int_ring_t* K, lp_integer_t* c, long x) {
  mpz_init_set_si(c, x);
  integer_ring_normalize(K, c);
}

static inline
void integer_construct_from_string(lp_int_ring_t* K, lp_integer_t* c, const char* x, int base) {
  mpz_init_set_str(c, x, base);
  integer_ring_normalize(K, c);
}

static inline
void integer_construct_copy(lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from) {
  mpz_init_set(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_assign(lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from) {
  mpz_set(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_assign_int(lp_int_ring_t* K, lp_integer_t* c, long from) {
  mpz_set_si(c, from);
  integer_ring_normalize(K, c);
}

static inline
void integer_destruct(lp_integer_t* c) {
  mpz_clear(c);
}

static inline
int integer_print(const lp_integer_t* c, FILE* out) {
  return mpz_out_str(out, 10, c);
}

static inline
int integer_print_matrix(const lp_integer_t* c, size_t m, size_t n, FILE* out) {
  size_t i, j;
  int len = 0;
  for (i = 0; i < m; ++ i) {
    for (j = 0; j < n; ++ j) {
      len += gmp_fprintf(out, "%4Zd", c + i*m + j);
    }
    len += fprintf(out, "\n");
  }
  return len;
}

static inline
char* integer_to_string(const lp_integer_t* c) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  integer_print(c, f);
  fclose(f);
  return str;
}


static inline
size_t integer_bits(const lp_integer_t* c) {
  return mpz_sizeinbase(c, 2);
}

static inline
long integer_to_int(const lp_integer_t* c) {
  return mpz_get_si(c);
}

static inline
int integer_is_prime(const lp_integer_t* c) {
  return mpz_probab_prime_p(c, 25);
}

static inline
int integer_is_zero(lp_int_ring_t* K, const lp_integer_t* c) {
  if (K) {
    lp_integer_t c_normalized;
    integer_construct_copy(K, &c_normalized, c);
    int sgn = mpz_sgn(&c_normalized);
    integer_destruct(&c_normalized);
    return sgn == 0;
  } else {
    return mpz_sgn(c) == 0;
  }
}

static inline
int integer_sgn(lp_int_ring_t* K, const lp_integer_t* c) {
  if (K) {
    lp_integer_t c_normalized;
    integer_construct_copy(K, &c_normalized, c);
    int sgn = mpz_sgn(&c_normalized);
    integer_destruct(&c_normalized);
    return sgn;
  } else {
    return mpz_sgn(c);
  }
}

static inline
int integer_cmp(lp_int_ring_t* K, const lp_integer_t* c, const lp_integer_t* to) {
  if (K) {
    lp_integer_t c_normalized, to_normalized;
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
int integer_cmp_int(lp_int_ring_t* K, const lp_integer_t* c, long to) {
  if (K) {
    lp_integer_t c_normalized, to_normalized;
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
int integer_divides(lp_int_ring_t* K, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  if (K) {
    // In a prime ring, it's always divisible
    if (K->is_prime) return integer_sgn(lp_Z, a);
    // Otherwise compute the gcd
    lp_integer_t gcd;
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
void integer_swap(lp_integer_t* a, lp_integer_t* b) {
  mpz_swap(a, b);
}

static inline
void integer_inc(lp_int_ring_t* K, lp_integer_t* a) {
  assert(integer_in_ring(K, a));
  lp_integer_t tmp;
  mpz_init(&tmp);
  mpz_add_ui(&tmp, a, 1);
  mpz_swap(&tmp, a);
  mpz_clear(&tmp);
  integer_ring_normalize(K, a);
}

static inline
void integer_dec(lp_int_ring_t* K, lp_integer_t* a) {
  assert(integer_in_ring(K, a));
  lp_integer_t tmp;
  mpz_init(&tmp);
  mpz_sub_ui(&tmp, a, 1);
  mpz_swap(&tmp, a);
  mpz_clear(&tmp);
  integer_ring_normalize(K, a);
}

static inline
void integer_add(lp_int_ring_t* K, lp_integer_t* sum, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_add(sum, a, b);
  integer_ring_normalize(K, sum);
}

static inline
void integer_sub(lp_int_ring_t* K, lp_integer_t* sub, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_sub(sub, a, b);
  integer_ring_normalize(K, sub);
}

static inline
void integer_neg(lp_int_ring_t* K, lp_integer_t* neg, const lp_integer_t* a) {
  assert(integer_in_ring(K, a));
  mpz_neg(neg, a);
  integer_ring_normalize(K, neg);
}

static inline
void integer_abs(lp_int_ring_t* K, lp_integer_t* abs, const lp_integer_t* a) {
  assert(integer_in_ring(K, a));
  mpz_abs(abs, a);
  integer_ring_normalize(K, abs);
}

static inline
void integer_inv(lp_int_ring_t* K, lp_integer_t* inv, const lp_integer_t* a) {
  assert(K);
  assert(integer_in_ring(K, a));
  int result = mpz_invert(inv, a, &K->M);
  assert(result);
  __unused(result);
  integer_ring_normalize(K, inv);
}

static inline
void integer_mul(lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_mul(product, a, b);
  integer_ring_normalize(K, product);
}

static inline
void integer_mul_int(lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, long b) {
  assert(integer_in_ring(K, a));
  mpz_mul_si(product, a, b);
  integer_ring_normalize(K, product);
}

static inline
void integer_mul_pow2(lp_int_ring_t* K, lp_integer_t* power, const lp_integer_t* a, unsigned n) {
  assert(integer_in_ring(K, a));
  assert(n > 0);
  mpz_mul_2exp(power, a, n);
  integer_ring_normalize(K, power);
}

static inline
void integer_pow(lp_int_ring_t* K, lp_integer_t* power, const lp_integer_t*a, unsigned n) {
  assert(integer_in_ring(K, a));
  if (K) {
    mpz_powm_ui(power, a, n, &K->M);
    integer_ring_normalize(K, power);
  } else {
    mpz_pow_ui(power, a, n);
  }
}

static inline
void integer_sqrt_Z(lp_integer_t* sqrt, const lp_integer_t* a) {
  mpz_sqrt(sqrt, a);
}

static inline
void integer_add_mul(lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, sum_product) && integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_addmul(sum_product, a, b);
  integer_ring_normalize(K, sum_product);
}

static inline
void integer_sub_mul(lp_int_ring_t* K, lp_integer_t* sub_product, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, sub_product) && integer_in_ring(K, a) && integer_in_ring(K, b));
  mpz_submul(sub_product, a, b);
  integer_ring_normalize(K, sub_product);
}

static inline
void integer_add_mul_int(lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, int b) {
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
void integer_div_exact(lp_int_ring_t* K, lp_integer_t* div, const lp_integer_t* a, const lp_integer_t* b) {
  assert(integer_in_ring(K, a) && integer_in_ring(K, b));
  if (K) {
    // Solving a = div*b (mod M). Let d = gcd(b, M) with extended gcd, we have
    // that c1*b+c2*M = d. Since d should divide a, the we get the solution
    // multiplying by a/d, obtaining (c1*a/d)*b = a.
    lp_integer_t c1, c2, gcd;
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
void integer_div_Z(lp_integer_t* div, const lp_integer_t* a, const lp_integer_t* b) {
  mpz_tdiv_q(div, a, b);
}

static inline
void integer_rem_Z(lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b) {
  mpz_tdiv_r(rem, a, b);
}

static inline
void integer_div_rem_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b) {
  mpz_tdiv_qr(div, rem, a, b);
}

static inline
void integer_div_rem_pow2_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, unsigned n) {
  mpz_tdiv_q_2exp(div, a, n);
  mpz_tdiv_r_2exp(rem, a, n);
}

static inline
void integer_gcd_Z(lp_integer_t* gcd, const lp_integer_t* a, const lp_integer_t* b) {
  mpz_gcd(gcd, a, b);
}

static inline
void integer_lcm_Z(lp_integer_t* lcm, const lp_integer_t* a, const lp_integer_t* b) {
  mpz_lcm(lcm, a, b);
}
