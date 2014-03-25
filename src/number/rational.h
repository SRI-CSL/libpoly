/*
 * rational_internal.h
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#pragma once

#include <rational.h>
#include <assert.h>

static inline
void rational_construct(rational_t* q) {
  mpq_init(q);
}

static inline
void rational_construct_from_int(rational_t* q, long a, unsigned long b) {
  mpq_init(q);
  mpq_set_si(q, a, b);
  mpq_canonicalize(q);
}

static inline
void rational_construct_from_integer(rational_t* q, const integer_t* from) {
  mpq_init(q);
  mpq_set_z(q, from);
}

static inline
void rational_construct_from_double(rational_t* q, double from) {
  mpq_init(q);
  mpq_set_d(q, from);
}

static inline
void rational_construct_from_dyadic(rational_t* q, const dyadic_rational_t* qd) {
  mpq_init(q);
  mpq_set_z(q, &qd->a);
  if (qd->n > 0) {
    mpq_div_2exp(q, q, qd->n);
  }
}

static inline
void rational_construct_copy(rational_t* q, const rational_t* from) {
  mpq_init(q);
  mpq_set(q, from);
}

static inline
void rational_assign(rational_t* q, const rational_t* from) {
  mpq_set(q, from);
}

static inline
void rational_assign_int(rational_t* q, long a, unsigned long b) {
  mpq_set_si(q, a, b);
  mpq_canonicalize(q);
}

static inline
void rational_destruct(rational_t* q) {
  mpq_clear(q);
}

static inline
int rational_print(const rational_t* c, FILE* out) {
  return mpq_out_str(out, 10, c);
}

static inline
char* rational_to_string(const rational_t* q) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  rational_print(q, f);
  fclose(f);
  return str;
}

static inline
double rational_to_double(const rational_t* q) {
  return mpq_get_d(q);
}

static inline
int rational_sgn(const rational_t* q) {
  return mpq_sgn(q);
}

static inline
int rational_cmp(const rational_t* q1, const rational_t* q2) {
  return mpq_cmp(q1, q2);
}

static inline
void rational_swap(rational_t* q1, rational_t* q2) {
  mpq_swap(q1, q2);
}

static inline
void rational_add(rational_t* sum, const rational_t* a, const rational_t* b) {
  mpq_add(sum, a, b);
}

static inline
void rational_add_integer(rational_t* sum, const rational_t* a, const integer_t* b) {
  rational_t b_rat;
  rational_construct_from_integer(&b_rat, b);
  mpq_add(sum, a, &b_rat);
  rational_destruct(&b_rat);
}

static inline
void rational_sub(rational_t* sub, const rational_t* a, const rational_t* b) {
  mpq_sub(sub, a, b);
}

static inline
void rational_neg(rational_t* neg, const rational_t* a) {
  mpq_neg(neg, a);
}

static inline
void rational_inv(rational_t* inv, const rational_t* a) {
  mpq_inv(inv, a);
}

static inline
void rational_mul(rational_t* mul, const rational_t* a, const rational_t* b) {
  mpq_mul(mul, a, b);
}

static inline
void rational_mul_2exp(rational_t* mul, const rational_t* a, unsigned n) {
  mpq_mul_2exp(mul, a, n);
}

static inline
void rational_pow(rational_t* pow, const rational_t*a, unsigned n) {
  rational_t result, tmp;
  rational_construct_from_int(&result, 1, 1);
  rational_construct_copy(&tmp, a);
  while (n) {
    if (n & 1) {
      rational_mul(&result, &result, &tmp);
    }
    rational_mul(&tmp, &tmp, &tmp);
    n >>= 1;
  }
  rational_swap(&result, pow);
  rational_destruct(&tmp);
  rational_destruct(&result);
}

static inline
void rational_div(rational_t* div, const rational_t* a, const rational_t* b) {
  mpq_div(div, a, b);
}

static inline
void rational_div_2exp(rational_t* div, const rational_t* a, unsigned n) {
  mpq_div_2exp(div, a, n);
}
