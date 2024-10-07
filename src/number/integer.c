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

#include "number/integer.h"

lp_int_ring_t* lp_int_ring_create(const lp_integer_t* M, int is_prime) {

  assert(mpz_sgn(M) > 0);

  lp_integer_t tmp;
  mpz_init(&tmp);

  lp_int_ring_t* K = malloc(sizeof(lp_int_ring_t));

  K->is_prime = is_prime;
  assert(!!is_prime == !!mpz_probab_prime_p(M, 25));
  K->ref_count = 1;

  mpz_init_set(&K->M, M);
  mpz_init(&K->ub);
  mpz_tdiv_q_2exp(&K->ub, M, 1);

  mpz_init(&K->lb);
  mpz_sub_ui(&K->lb, M, 1); // M-1
  mpz_tdiv_q_2exp(&tmp, &K->lb, 1); // [M-1/2]
  mpz_neg(&K->lb, &tmp);

  mpz_clear(&tmp);

  return K;
}

void lp_int_ring_destroy(lp_int_ring_t* K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  mpz_clear(&nonconst->M);
  mpz_clear(&nonconst->ub);
  mpz_clear(&nonconst->lb);
  free(nonconst);
}

void lp_int_ring_attach(lp_int_ring_t* K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    ++ nonconst->ref_count;
  }
}

void lp_int_ring_detach(lp_int_ring_t* K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    assert(nonconst->ref_count > 0);
    if (--nonconst->ref_count == 0) {
      lp_int_ring_destroy(nonconst);
    }
  }
}

int lp_int_ring_equal(const lp_int_ring_t* K1, const lp_int_ring_t* K2) {
  if (K1 == K2) {
    return 1;
  } else if (K1 && K2) {
    return mpz_cmp(&K1->M, &K2->M) == 0;
  } else {
    return 0;
  }
}


int lp_int_ring_print(const lp_int_ring_t* K, FILE* out) {
  int len = 0;
  len += fprintf(out, "Z");
  if (K) {
    len += fprintf(out, " mod ");
    len += integer_print(&K->M, out);
  }
  return len;
}

char* lp_int_ring_to_string(const lp_int_ring_t* K) {
  struct u_memstream mem;
  char* str = 0;
  size_t size = 0;
  u_memstream_open(&mem, &str, &size);
  FILE* f = u_memstream_get(&mem);
  lp_int_ring_print(K, f);
  u_memstream_close(&mem);
  return str;
}

lp_int_ring_t* lp_Z = 0;

void lp_integer_construct(lp_integer_t* c) {
  integer_construct(c);
}

void lp_integer_construct_from_rational(const lp_int_ring_t* K, lp_integer_t* c, const lp_rational_t* q) {
  integer_construct_from_rational(K, c, q);
}

void lp_integer_construct_from_int(const lp_int_ring_t* K, lp_integer_t* c, long x) {
  integer_construct_from_int(K, c, x);
}

void lp_integer_construct_from_string(const lp_int_ring_t* K, lp_integer_t* c, const char* x, int base) {
  integer_construct_from_string(K, c, x, base);
}

void lp_integer_construct_copy(const lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from) {
  integer_construct_copy(K, c, from);
}

void lp_integer_assign(const lp_int_ring_t* K, lp_integer_t* c, const lp_integer_t* from) {
  integer_assign(K, c, from);
}

void lp_integer_assign_int(const lp_int_ring_t* K, lp_integer_t* c, long x) {
  integer_assign_int(K, c, x);
}

void lp_integer_destruct(lp_integer_t* c) {
  integer_destruct(c);
}

int lp_integer_print(const lp_integer_t* c, FILE* out) {
  return integer_print(c, out);
}

size_t lp_integer_bits(const lp_integer_t* c) {
  return integer_bits(c);
}

int lp_integer_print_matrix(const lp_integer_t* c, size_t m, size_t n, FILE* out) {
  return integer_print_matrix(c, m, n, out);
}

char* lp_integer_to_string(const lp_integer_t* c) {
  return integer_to_string(c);
}

int lp_integer_fits_int(const lp_integer_t* c) {
  return integer_fits_int(c);
}

long lp_integer_to_int(const lp_integer_t* c) {
  return integer_to_int(c);
}

double lp_integer_to_double(const lp_integer_t* c) {
  return integer_to_double(c);
}

int lp_integer_is_prime(const lp_integer_t* c) {
  return integer_is_prime(c);
}

int lp_integer_is_zero(const lp_int_ring_t* K, const lp_integer_t* c) {
  return integer_is_zero(K, c);
}

int lp_integer_in_ring(const lp_int_ring_t* K, const lp_integer_t* c) {
  return integer_in_ring(K, c);
}

int lp_integer_sgn(const lp_int_ring_t* K, const lp_integer_t* c) {
  return integer_sgn(K, c);
}

int lp_integer_cmp(const lp_int_ring_t* K, const lp_integer_t* c, const lp_integer_t* to) {
  return integer_cmp(K, c, to);
}

int lp_integer_cmp_int(const lp_int_ring_t* K, const lp_integer_t* c, long to) {
  return integer_cmp_int(K, c, to);
}

int lp_integer_divides(const lp_int_ring_t* K, const lp_integer_t* a, const lp_integer_t* b) {
  return integer_divides(K, a, b);
}

void lp_integer_swap(lp_integer_t* a, lp_integer_t* b) {
  integer_swap(a, b);
}

void lp_integer_inc(const lp_int_ring_t* K, lp_integer_t* a) {
  integer_inc(K, a);
}

void lp_integer_dec(const lp_int_ring_t* K, lp_integer_t* a) {
  integer_dec(K, a);
}

void lp_integer_add(const lp_int_ring_t* K, lp_integer_t* sum, const lp_integer_t* a, const lp_integer_t* b) {
  integer_add(K, sum, a, b);
}

void lp_integer_sub(const lp_int_ring_t* K, lp_integer_t* sub, const lp_integer_t* a, const lp_integer_t* b) {
  integer_sub(K, sub, a, b);
}

void lp_integer_neg(const lp_int_ring_t* K, lp_integer_t* neg, const lp_integer_t* a) {
  integer_neg(K, neg, a);
}

void lp_integer_abs(const lp_int_ring_t* K, lp_integer_t* abs, const lp_integer_t* a) {
  integer_abs(K, abs, a);
}

void lp_integer_inv(const lp_int_ring_t* K, lp_integer_t* inv, const lp_integer_t* a) {
  integer_inv(K, inv, a);
}

void lp_integer_mul(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, const lp_integer_t* b) {
  integer_mul(K, product, a, b);
}

void lp_integer_mul_int(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, long b) {
  integer_mul_int(K, product, a, b);
}

void lp_integer_mul_pow2(const lp_int_ring_t* K, lp_integer_t* product, const lp_integer_t* a, unsigned n) {
  integer_mul_pow2(K, product, a, n);
}

void lp_integer_pow(const lp_int_ring_t* K, lp_integer_t* pow, const lp_integer_t *a, unsigned n) {
  integer_pow(K, pow, a, n);
}

void lp_integer_sqrt_Z(lp_integer_t* sqrt, const lp_integer_t* a) {
  integer_sqrt_Z(sqrt, a);
}

void lp_integer_add_mul(const lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, const lp_integer_t* b) {
  integer_add_mul(K, sum_product, a, b);
}

void lp_integer_add_mul_int(const lp_int_ring_t* K, lp_integer_t* sum_product, const lp_integer_t* a, int b) {
  integer_add_mul_int(K, sum_product, a, b);
}

void lp_integer_sub_mul(const lp_int_ring_t* K, lp_integer_t* sub_product, const lp_integer_t* a, const lp_integer_t* b) {
  integer_sub_mul(K, sub_product, a, b);
}

void lp_integer_div_exact(const lp_int_ring_t* K, lp_integer_t* div_Z, const lp_integer_t* a, const lp_integer_t* b) {
  integer_div_exact(K, div_Z, a, b);
}

void lp_integer_div_Z(lp_integer_t* div, const lp_integer_t* a, const lp_integer_t* b) {
  integer_div_Z(div, a, b);
}

void lp_integer_rem_Z(lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b) {
  integer_rem_Z(rem, a, b);
}

void lp_integer_div_rem_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, const lp_integer_t* b) {
  integer_div_rem_Z(div, rem, a, b);
}

void lp_integer_div_rem_pow2_Z(lp_integer_t* div, lp_integer_t* rem, const lp_integer_t* a, unsigned n) {
  integer_div_rem_pow2_Z(div, rem, a, n);
}

void lp_integer_gcd_Z(lp_integer_t* gcd, const lp_integer_t* a, const lp_integer_t* b) {
  integer_gcd_Z(gcd, a, b);
}

void lp_integer_lcm_Z(lp_integer_t* lcm, const lp_integer_t* a, const lp_integer_t* b) {
  integer_lcm_Z(lcm, a, b);
}

size_t lp_integer_hash(const lp_integer_t* a) {
  return integer_hash(a);
}
