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

#include "number/rational.h"
#include "number/integer.h"

#include <assert.h>

/** Normalize q, so that a and 2^n have no common factors. */
static inline
void dyadic_rational_normalize(lp_dyadic_rational_t* q) {
  if (mpz_sgn(&q->a) == 0) {
    q->n = 0;
  } else if (q->n > 0) {
    unsigned long first1 = mpz_scan1(&q->a, 0);
    if (first1 > 0) {
      unsigned long min = first1 >= q->n ? q->n : first1;
      q->n -= min;
      mpz_div_2exp(&q->a, &q->a, min);
    }
  }
}

/** Check if normalized as above */
static inline
int dyadic_rational_is_normalized(const lp_dyadic_rational_t* q) {
  return (mpz_sgn(&q->a) == 0 && q->n == 0) || (mpz_scan1(&q->a, 0) == 0 || q->n == 0);
}

static inline
void dyadic_rational_construct(lp_dyadic_rational_t* q) {
  mpz_init(&q->a);
  q->n = 0;
}

static inline
void dyadic_rational_construct_from_int(lp_dyadic_rational_t* q, long a, unsigned long n) {
  mpz_init_set_si(&q->a, a);
  q->n = n;
  dyadic_rational_normalize(q);
}

static inline
void dyadic_rational_construct_from_double(lp_dyadic_rational_t* q, double x) {
  lp_rational_t tmp;
  mpq_init(&tmp);
  mpq_set_d(&tmp, x);
  mpz_init_set(&q->a, &tmp._mp_num);
  q->n = mpz_scan1(&tmp._mp_den, 0);
  mpq_clear(&tmp);
  assert(dyadic_rational_is_normalized(q));
}

static inline
void dyadic_rational_construct_from_integer(lp_dyadic_rational_t* q, const lp_integer_t* from) {
  mpz_init_set(&q->a, from);
  q->n = 0;
  dyadic_rational_normalize(q);
}

static inline
void dyadic_rational_construct_copy(lp_dyadic_rational_t* q, const lp_dyadic_rational_t* from) {
  assert(dyadic_rational_is_normalized(from));
  mpz_init_set(&q->a, &from->a);
  q->n = from->n;
}

static inline
void dyadic_rational_assign(lp_dyadic_rational_t* q, const lp_dyadic_rational_t* from) {
  assert(dyadic_rational_is_normalized(from));
  mpz_set(&q->a, &from->a);
  q->n = from->n;
}

static inline
void dyadic_rational_assign_int(lp_dyadic_rational_t* q, long a, unsigned long n) {
  mpz_set_si(&q->a, a);
  q->n = n;
  dyadic_rational_normalize(q);
}

static inline
void dyadic_rational_destruct(lp_dyadic_rational_t* q) {
  mpz_clear(&q->a);
}

static inline
int dyadic_rational_print(const lp_dyadic_rational_t* dq, FILE* out) {
  lp_rational_t q;
  rational_construct_from_dyadic(&q, dq);
  int ret = rational_print(&q, out);
  rational_destruct(&q);
  return ret;
}

static inline
char* dyadic_rational_to_string(const lp_dyadic_rational_t* q) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  dyadic_rational_print(q, f);
  fclose(f);
  return str;
}

static inline
double dyadic_rational_to_double(const lp_dyadic_rational_t* dq) {
  lp_rational_t q;
  rational_construct_from_dyadic(&q, dq);
  double ret = rational_to_double(&q);
  rational_destruct(&q);
  return ret;
}

static inline
int dyadic_rational_sgn(const lp_dyadic_rational_t* q) {
  assert(dyadic_rational_is_normalized(q));
  return mpz_sgn(&q->a);
}

static inline
int dyadic_rational_cmp(const lp_dyadic_rational_t* q1, const lp_dyadic_rational_t* q2) {
  assert(dyadic_rational_is_normalized(q1));
  assert(dyadic_rational_is_normalized(q2));
  int q1_sgn = dyadic_rational_sgn(q1);
  int q2_sgn = dyadic_rational_sgn(q2);
  if (q1_sgn == q2_sgn) {
    if (q1_sgn == 0) return 0;
    if (q1->n == q2->n) {
      // sgn(a - b)
      return mpz_cmp(&q1->a, &q2->a);
    } else if (q1->n > q2->n) {
      // sgn(a/2^n - b/2^m) = sgn(a - b*2^(n-m))
      lp_integer_t q2tmp;
      mpz_init(&q2tmp);
      mpz_mul_2exp(&q2tmp, &q2->a, q1->n - q2->n);
      int cmp = mpz_cmp(&q1->a, &q2tmp);
      mpz_clear(&q2tmp);
      return cmp;
    } else {
      // sgn(a/2^n - b/2^m) = sgn(a*2^(m-n) - b)
      lp_integer_t q1tmp;
      mpz_init(&q1tmp);
      mpz_mul_2exp(&q1tmp, &q1->a, q2->n - q1->n);
      int cmp = mpz_cmp(&q1tmp, &q2->a);
      mpz_clear(&q1tmp);
      return cmp;
    }
  } else {
    return q1_sgn - q2_sgn;
  }
}

static inline
int dyadic_rational_cmp_integer(const lp_dyadic_rational_t* q1, const lp_integer_t* q2) {
  lp_dyadic_rational_t q2_dy_rat;
  dyadic_rational_construct_from_integer(&q2_dy_rat, q2);
  int cmp = dyadic_rational_cmp(q1, &q2_dy_rat);
  dyadic_rational_destruct(&q2_dy_rat);
  return cmp;
}

static inline
void dyadic_rational_swap(lp_dyadic_rational_t* q1, lp_dyadic_rational_t* q2) {
  assert(dyadic_rational_is_normalized(q1));
  assert(dyadic_rational_is_normalized(q2));
  mpz_swap(&q1->a, &q2->a);
  long tmp = q1->n;
  q1->n = q2->n;
  q2->n = tmp;
}

static inline
void dyadic_rational_add(lp_dyadic_rational_t* sum, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  assert(dyadic_rational_is_normalized(a));
  assert(dyadic_rational_is_normalized(b));

  if (a->n == b->n) {
    mpz_add(&sum->a, &a->a, &b->a);
    sum->n = a->n;
  } else if (a->n > b->n) {
    // m > n
    // a/2^m + b/2^n = (a + 2^(m-n)*b)/2^m
    lp_integer_t b2n;
    mpz_init(&b2n);
    mpz_mul_2exp(&b2n, &b->a, a->n - b->n);
    mpz_add(&sum->a, &a->a, &b2n);
    mpz_clear(&b2n);
    sum->n = a->n;
  } else {
    // m < n
    // a/2^m + b/2^n = (a*2^(n-m) + b)/2^n
    lp_integer_t a2n;
    mpz_init(&a2n);
    mpz_mul_2exp(&a2n, &a->a, b->n - a->n);
    mpz_add(&sum->a, &a2n, &b->a);
    mpz_clear(&a2n);
    sum->n = b->n;
  }

  dyadic_rational_normalize(sum);
}

static inline
void dyadic_rational_add_integer(lp_dyadic_rational_t* sum, const lp_dyadic_rational_t* a, const lp_integer_t* b) {
  assert(dyadic_rational_is_normalized(a));
  if (a->n > 0) {
    // a/2^n + b = (a + b*2^n)/2^n
    lp_integer_t b2n;
    mpz_init(&b2n);
    mpz_mul_2exp(&b2n, b, a->n);
    mpz_add(&sum->a, &a->a, &b2n);
    mpz_clear(&b2n);
  } else {
    mpz_add(&sum->a, &a->a, b);
  }
  sum->n = a->n;
  dyadic_rational_normalize(sum);
}

static inline
void dyadic_rational_sub(lp_dyadic_rational_t* sub, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  assert(dyadic_rational_is_normalized(a));
  assert(dyadic_rational_is_normalized(b));

  if (a->n == b->n) {
    mpz_sub(&sub->a, &a->a, &b->a);
    sub->n = a->n;
  } else if (a->n > b->n) {
    // m > n
    // a/2^m - b/2^n = (a - 2^(m-n)*b)/2^m
    lp_integer_t b2n;
    mpz_init(&b2n);
    mpz_mul_2exp(&b2n, &b->a, a->n - b->n);
    mpz_sub(&sub->a, &a->a, &b2n);
    mpz_clear(&b2n);
    sub->n = a->n;
  } else {
    // m < n
    // a/2^m - b/2^n = (a^(n-m) - b)/2^n
    lp_integer_t a2n;
    mpz_init(&a2n);
    mpz_mul_2exp(&a2n, &a->a, b->n - a->n);
    mpz_sub(&sub->a, &a2n, &b->a);
    mpz_clear(&a2n);
    sub->n = b->n;
  }

  dyadic_rational_normalize(sub);
}

static inline
void dyadic_rational_neg(lp_dyadic_rational_t* neg, const lp_dyadic_rational_t* a) {
  assert(dyadic_rational_is_normalized(a));
  mpz_neg(&neg->a, &a->a);
}

static inline
void dyadic_rational_mul(lp_dyadic_rational_t* mul, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  assert(dyadic_rational_is_normalized(a));
  assert(dyadic_rational_is_normalized(b));
  mpz_mul(&mul->a, &a->a, &b->a);
  mul->n = a->n + b->n;
  dyadic_rational_normalize(mul);
}

static inline
void dyadic_rational_mul_2exp(lp_dyadic_rational_t* mul, const lp_dyadic_rational_t* a, unsigned long n) {
  assert(dyadic_rational_is_normalized(a));
  mpz_set(&mul->a, &a->a);
  // a/2^n + 2^n
  if (a->n >= n) {
    mul->n = a->n - n;
  } else {
    mpz_mul_2exp(&mul->a, &a->a, n - mul->n);
    mul->n = 0;
  }
}

static inline
void dyadic_rational_pow(lp_dyadic_rational_t* pow, const lp_dyadic_rational_t* a, unsigned long n) {
  assert(dyadic_rational_is_normalized(a));
  mpz_pow_ui(&pow->a, &a->a, n);
  pow->n = a->n * n;
}

static inline
void dyadic_rational_div_2exp(lp_dyadic_rational_t* div, const lp_dyadic_rational_t* a, unsigned long n) {
  assert(dyadic_rational_is_normalized(a));
  mpz_set(&div->a, &a->a);
  div->n = a->n + n;
  dyadic_rational_normalize(div);
}

static inline
void dyadic_rational_get_num(const lp_dyadic_rational_t* q, lp_integer_t* num) {
  integer_assign(lp_Z, num, &q->a);
}

static inline
void dyadic_rational_get_den(const lp_dyadic_rational_t* q, lp_integer_t* den) {
  integer_assign_int(lp_Z, den, 1);
  integer_mul_pow2(lp_Z, den, den, q->n);
}

static inline
int dyadic_rational_get_distance_size(const lp_dyadic_rational_t* lower, const lp_dyadic_rational_t* upper) {

  assert(dyadic_rational_cmp(lower, upper) < 0);

  int size;

  if (lower->n == upper->n) {
    // size([l/2^n, u/2^n]) = size([l,u]) - n
    lp_integer_t diff;
    integer_construct(&diff);
    integer_sub(lp_Z, &diff, &upper->a, &lower->a);
    size = integer_log2_abs(&diff);
    size -= lower->n;
    integer_destruct(&diff);
  } else if (lower->n > upper->n) {
    // n1 > n2
    // size([l/2^n1, u/2^n2]) = log2( u*2^(n1 - n2) - l) / 2^n1)
    lp_integer_t diff;
    integer_construct(&diff);
    integer_mul_pow2(lp_Z, &diff, &upper->a, lower->n - upper->n);
    integer_sub(lp_Z, &diff, &diff, &lower->a);
    size = integer_log2_abs(&diff);
    size -= lower->n;
    integer_destruct(&diff);
  } else {
    // n1 < n2
    // size([l/2^n1, u/2^n2]) = log2( u - l*2^(n2 - n1)) / 2^n2)
    lp_integer_t diff;
    integer_construct(&diff);
    integer_mul_pow2(lp_Z, &diff, &lower->a, upper->n - lower ->n);
    integer_sub(lp_Z, &diff, &upper->a, &diff);
    size = integer_log2_abs(&diff);
    size -= upper->n;
    integer_destruct(&diff);
  }

  return size;
}

static inline
void dyadic_rational_ceiling(const lp_dyadic_rational_t* a, lp_dyadic_rational_t* ceil) {
  // Just divide a with 2^n
  if (a->n > 0) {
    integer_div_ceiling_pow2(&ceil->a, &a->a, a->n);
    ceil->n = 0;
  } else {
    dyadic_rational_assign(ceil, a);
  }
}

static inline
void dyadic_rational_floor(const lp_dyadic_rational_t* a, lp_dyadic_rational_t* floor) {
  // Just divide a with 2^n
  if (a->n > 0) {
    integer_div_floor_pow2(&floor->a, &a->a, a->n);
    floor->n = 0;
  } else {
    dyadic_rational_assign(floor, a);
  }
}

static inline
void dyadic_rational_ceiling_int(const lp_dyadic_rational_t* a, lp_integer_t* ceil) {
  // Just divide a with 2^n
  if (a->n > 0) {
    integer_div_ceiling_pow2(ceil, &a->a, a->n);
  } else {
    integer_assign(lp_Z, ceil, &a->a);
  }
}

static inline
void dyadic_rational_floor_int(const lp_dyadic_rational_t* a, lp_integer_t* floor) {
  // Just divide a with 2^n
  if (a->n > 0) {
    integer_div_floor_pow2(floor, &a->a, a->n);
  } else {
    integer_assign(lp_Z, floor, &a->a);
  }
}

static inline
int dyadic_rational_is_integer(const lp_dyadic_rational_t* a) {
  return a->n == 0;
}

/** Returns a dyadic rational in the interval (a, b) */
void dyadic_rational_get_value_between(lp_dyadic_rational_t* v, const lp_rational_t* a, const lp_rational_t* b);

