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
void rational_construct(lp_rational_t* q) {
  mpq_init(q);
}

static inline
void rational_construct_from_int(lp_rational_t* q, long a, unsigned long b) {
  mpq_init(q);
  mpq_set_si(q, a, b);
  mpq_canonicalize(q);
}

static inline
void rational_construct_from_integer(lp_rational_t* q, const lp_integer_t* from) {
  mpq_init(q);
  mpq_set_z(q, from);
}

static inline
void rational_construct_from_double(lp_rational_t* q, double from) {
  mpq_init(q);
  mpq_set_d(q, from);
}

static inline
void rational_construct_from_dyadic(lp_rational_t* q, const lp_dyadic_rational_t* qd) {
  mpq_init(q);
  mpq_set_z(q, &qd->a);
  if (qd->n > 0) {
    mpq_div_2exp(q, q, qd->n);
  }
}

static inline
void rational_construct_copy(lp_rational_t* q, const lp_rational_t* from) {
  mpq_init(q);
  mpq_set(q, from);
}

static inline
void rational_assign(lp_rational_t* q, const lp_rational_t* from) {
  mpq_set(q, from);
}

static inline
void rational_assign_int(lp_rational_t* q, long a, unsigned long b) {
  mpq_set_si(q, a, b);
  mpq_canonicalize(q);
}

static inline
void rational_destruct(lp_rational_t* q) {
  mpq_clear(q);
}

static inline
int rational_print(const lp_rational_t* c, FILE* out) {
  return mpq_out_str(out, 10, c);
}

static inline
char* rational_to_string(const lp_rational_t* q) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  rational_print(q, f);
  fclose(f);
  return str;
}

static inline
double rational_to_double(const lp_rational_t* q) {
  return mpq_get_d(q);
}

static inline
int rational_sgn(const lp_rational_t* q) {
  return mpq_sgn(q);
}

static inline
int rational_cmp(const lp_rational_t* q1, const lp_rational_t* q2) {
  return mpq_cmp(q1, q2);
}

static inline
int rational_cmp_dyadic_rational(const lp_rational_t* q1, const lp_dyadic_rational_t* q2) {
  lp_rational_t q2_rat;
  rational_construct_from_dyadic(&q2_rat, q2);
  int cmp = rational_cmp(q1, &q2_rat);
  rational_destruct(&q2_rat);
  return cmp;
}

static inline
int rational_cmp_integer(const lp_rational_t* q1, const lp_integer_t* q2) {
  lp_rational_t q2_rat;
  rational_construct_from_integer(&q2_rat, q2);
  int cmp = rational_cmp(q1, &q2_rat);
  rational_destruct(&q2_rat);
  return cmp;
}

static inline
void rational_swap(lp_rational_t* q1, lp_rational_t* q2) {
  mpq_swap(q1, q2);
}

static inline
void rational_add(lp_rational_t* sum, const lp_rational_t* a, const lp_rational_t* b) {
  mpq_add(sum, a, b);
}

static inline
void rational_add_integer(lp_rational_t* sum, const lp_rational_t* a, const lp_integer_t* b) {
  lp_rational_t b_rat;
  rational_construct_from_integer(&b_rat, b);
  mpq_add(sum, a, &b_rat);
  rational_destruct(&b_rat);
}

static inline
void rational_sub(lp_rational_t* sub, const lp_rational_t* a, const lp_rational_t* b) {
  mpq_sub(sub, a, b);
}

static inline
void rational_neg(lp_rational_t* neg, const lp_rational_t* a) {
  mpq_neg(neg, a);
}

static inline
void rational_inv(lp_rational_t* inv, const lp_rational_t* a) {
  mpq_inv(inv, a);
}

static inline
void rational_mul(lp_rational_t* mul, const lp_rational_t* a, const lp_rational_t* b) {
  mpq_mul(mul, a, b);
}

static inline
void rational_mul_2exp(lp_rational_t* mul, const lp_rational_t* a, unsigned n) {
  mpq_mul_2exp(mul, a, n);
}

static inline
void rational_pow(lp_rational_t* pow, const lp_rational_t*a, unsigned n) {
  lp_rational_t result, tmp;
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
void rational_div(lp_rational_t* div, const lp_rational_t* a, const lp_rational_t* b) {
  mpq_div(div, a, b);
}

static inline
void rational_div_2exp(lp_rational_t* div, const lp_rational_t* a, unsigned n) {
  mpq_div_2exp(div, a, n);
}

static inline
void rational_get_num(const lp_rational_t* q, lp_integer_t* num) {
  mpq_get_num(num, q);
}

static inline
void rational_get_den(const lp_rational_t* q, lp_integer_t* den) {
  mpq_get_den(den, q);
}

