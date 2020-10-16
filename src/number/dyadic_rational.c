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

#include "number/dyadic_rational.h"
#include "utils/hash.h"

void lp_dyadic_rational_construct(lp_dyadic_rational_t* q) {
  dyadic_rational_construct(q);
}

void lp_dyadic_rational_construct_from_int(lp_dyadic_rational_t* q, long a, unsigned long n) {
  dyadic_rational_construct_from_int(q, a, n);
}


void lp_dyadic_rational_construct_from_integer(lp_dyadic_rational_t* q, const lp_integer_t* a) {
  dyadic_rational_construct_from_integer(q, a);
}

void lp_dyadic_rational_construct_from_double(lp_dyadic_rational_t* q, double x) {
  dyadic_rational_construct_from_double(q, x);
}

void lp_dyadic_rational_construct_copy(lp_dyadic_rational_t* q, const lp_dyadic_rational_t* from) {
  dyadic_rational_construct_copy(q, from);
}


void lp_dyadic_rational_assign(lp_dyadic_rational_t* q, const lp_dyadic_rational_t* from) {
  dyadic_rational_assign(q, from);
}


void lp_dyadic_rational_assign_int(lp_dyadic_rational_t* q, long a, unsigned long n) {
  dyadic_rational_assign_int(q, a, n);
}


void lp_dyadic_rational_destruct(lp_dyadic_rational_t* q) {
  dyadic_rational_destruct(q);
}


int lp_dyadic_rational_print(const lp_dyadic_rational_t* c, FILE* out) {
  return dyadic_rational_print(c, out);
}


char* lp_dyadic_rational_to_string(const lp_dyadic_rational_t* q) {
  return dyadic_rational_to_string(q);
}


double lp_dyadic_rational_to_double(const lp_dyadic_rational_t* q) {
  return dyadic_rational_to_double(q);
}


int lp_dyadic_rational_sgn(const lp_dyadic_rational_t* q) {
  return dyadic_rational_sgn(q);
}


int lp_dyadic_rational_cmp(const lp_dyadic_rational_t* q1, const lp_dyadic_rational_t* q2) {
  return dyadic_rational_cmp(q1, q2);
}

int lp_dyadic_rational_cmp_integer(const lp_dyadic_rational_t* q1, const lp_integer_t* q2) {
  return dyadic_rational_cmp_integer(q1, q2);
}

int lp_dyadic_rational_cmp_rational(const lp_dyadic_rational_t* q1, const lp_rational_t* q2) {
  return -rational_cmp_dyadic_rational(q2, q1);
}

void lp_dyadic_rational_swap(lp_dyadic_rational_t* q1, lp_dyadic_rational_t* q2) {
  dyadic_rational_swap(q1, q2);
}

void lp_dyadic_rational_add(lp_dyadic_rational_t* sum, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  dyadic_rational_add(sum, a, b);
}


void lp_dyadic_rational_add_integer(lp_dyadic_rational_t* sum, const lp_dyadic_rational_t* a, const lp_integer_t* b) {
  dyadic_rational_add_integer(sum, a, b);
}


void lp_dyadic_rational_sub(lp_dyadic_rational_t* sub, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  dyadic_rational_sub(sub, a, b);
}


void lp_dyadic_rational_neg(lp_dyadic_rational_t* neg, const lp_dyadic_rational_t* a) {
  dyadic_rational_neg(neg, a);
}


void lp_dyadic_rational_mul(lp_dyadic_rational_t* mul, const lp_dyadic_rational_t* a, const lp_dyadic_rational_t* b) {
  dyadic_rational_mul(mul, a, b);
}


void lp_dyadic_rational_mul_2exp(lp_dyadic_rational_t* mul, const lp_dyadic_rational_t* a, unsigned long n) {
  dyadic_rational_mul_2exp(mul, a, n);
}


void lp_dyadic_rational_pow(lp_dyadic_rational_t* pow, const lp_dyadic_rational_t* a, unsigned long n) {
  dyadic_rational_pow(pow, a, n);
}

void lp_dyadic_rational_div_2exp(lp_dyadic_rational_t* div, const lp_dyadic_rational_t* a, unsigned long n) {
  dyadic_rational_div_2exp(div, a, n);
}

void lp_dyadic_rational_get_num(const lp_dyadic_rational_t* q, lp_integer_t* num) {
  dyadic_rational_get_num(q, num);
}

void lp_dyadic_rational_get_den(const lp_dyadic_rational_t* q, lp_integer_t* den) {
  dyadic_rational_get_den(q, den);
}

int lp_dyadic_rational_is_integer(const lp_dyadic_rational_t* q) {
  return dyadic_rational_is_integer(q);
}

void lp_dyadic_rational_ceiling(const lp_dyadic_rational_t* q, lp_integer_t* q_ceiling) {
  dyadic_rational_ceiling_int(q, q_ceiling);
}

void lp_dyadic_rational_floor(const lp_dyadic_rational_t* q, lp_integer_t* q_floor) {
  dyadic_rational_floor_int(q, q_floor);
}

void dyadic_rational_get_value_between(lp_dyadic_rational_t* v, const lp_rational_t* a, const lp_rational_t* b) {

  int cmp;
  lp_dyadic_rational_t result;

  assert(rational_cmp(a, b) < 0);

  // a <= a_ub < b_lb <= b
  // m = (a_ub + b_lb)/2 so a < m < b
  lp_rational_t m_q;
  rational_construct(&m_q);
  rational_add(&m_q, a, b);
  rational_div_2exp(&m_q, &m_q, 1);

  // floor(m) <= m <= ceil(m)
  lp_integer_t m_floor, m_ceil;
  integer_construct(&m_floor);
  rational_floor(&m_q, &m_floor);
  integer_construct_copy(lp_Z, &m_ceil, &m_floor);
  integer_inc(lp_Z, &m_ceil);

  // If a < m_floor (or equal and bound allows it), we can take this value
  cmp = rational_cmp_integer(a, &m_floor);
  if (cmp < 0) {
    lp_dyadic_rational_construct_from_integer(&result, &m_floor);
  } else {
    // If m_ceil < b_lb (or equal and bound allows it), we can take this value
    cmp = rational_cmp_integer(b, &m_ceil);
    if (cmp > 0) {
      lp_dyadic_rational_construct_from_integer(&result, &m_ceil);
    } else {

      lp_dyadic_rational_t m_dy, lb, ub;
      lp_dyadic_rational_construct(&m_dy);
      lp_dyadic_rational_construct_from_integer(&lb, &m_floor);
      lp_dyadic_rational_construct_from_integer(&ub, &m_ceil);

      // We have to do the search
      for (;;) {

        // always: lb < a < b < ub
        dyadic_rational_add(&m_dy, &lb, &ub);
        dyadic_rational_div_2exp(&m_dy, &m_dy, 1);

        // lb < m < a => move lb to m
        cmp = rational_cmp_dyadic_rational(a, &m_dy);
        if (cmp >= 0) {
          dyadic_rational_swap(&m_dy, &lb);
          continue;
        }

        // b < m < ub => move ub to m
        cmp = rational_cmp_dyadic_rational(b, &m_dy);
        if (cmp <= 0) {
          dyadic_rational_swap(&ub, &m_dy);
          continue;
        }

        // Got it l <= m <= u
        dyadic_rational_construct_copy(&result, &m_dy);
        break;
      }

      dyadic_rational_destruct(&lb);
      dyadic_rational_destruct(&ub);
      dyadic_rational_destruct(&m_dy);
    }
  }

  // Store result
  dyadic_rational_swap(v, &result);

  // Remove temps
  dyadic_rational_destruct(&result);
  integer_destruct(&m_ceil);
  integer_destruct(&m_floor);
  rational_destruct(&m_q);
}

size_t lp_dyadic_rational_hash(const lp_dyadic_rational_t* q) {
  size_t h1 = lp_integer_hash(&q->a);
  return hash_combine(h1, q->n);
}

