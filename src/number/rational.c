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

#include "number/rational.h"

void lp_rational_construct(lp_rational_t* q) {
  rational_construct(q);
}

void lp_rational_construct_from_div(lp_rational_t* q, const lp_integer_t* num, const lp_integer_t* den) {
  rational_construct_from_div(q, num, den);
}

void lp_rational_construct_from_int(lp_rational_t* q, long a, unsigned long b) {
  rational_construct_from_int(q, a, b);
}

void lp_rational_construct_from_integer(lp_rational_t* q, const lp_integer_t* a) {
  rational_construct_from_integer(q, a);
}

void lp_rational_construct_from_double(lp_rational_t* q, double x) {
  rational_construct_from_double(q, x);
}

void lp_rational_construct_from_dyadic(lp_rational_t* q, const lp_dyadic_rational_t* qd) {
  rational_construct_from_dyadic(q, qd);
}

void lp_rational_construct_copy(lp_rational_t* q, const lp_rational_t* from) {
  rational_construct_copy(q, from);
}

void lp_rational_assign(lp_rational_t* q, const lp_rational_t* from) {
  rational_assign(q, from);
}

void lp_rational_assign_int(lp_rational_t* q, long a, unsigned long b) {
  rational_assign_int(q, a, b);
}

void lp_rational_destruct(lp_rational_t* q) {
  rational_destruct(q);
}

int lp_rational_print(const lp_rational_t* c, FILE* out) {
  return rational_print(c, out);
}

#if _XOPEN_SOURCE >= 700 || _POSIX_C_SOURCE >= 200809L
char* lp_rational_to_string(const lp_rational_t* q) {
  return rational_to_string(q);
}
#endif


double lp_rational_to_double(const lp_rational_t* q) {
  return rational_to_double(q);
}

int lp_rational_sgn(const lp_rational_t* q) {
  return rational_sgn(q);
}

int lp_rational_cmp(const lp_rational_t* q1, const lp_rational_t* q2) {
  return rational_cmp(q1, q2);
}

int lp_rational_cmp_dyadic_rational(const lp_rational_t* q1, const lp_dyadic_rational_t* q2) {
  return rational_cmp_dyadic_rational(q1, q2);
}

int lp_rational_cmp_integer(const lp_rational_t* q1, const lp_integer_t* q2) {
  return rational_cmp_integer(q1, q2);
}


void lp_rational_swap(lp_rational_t* q1, lp_rational_t* q2) {
  rational_swap(q1, q2);
}

void lp_rational_add(lp_rational_t* sum, const lp_rational_t* a, const lp_rational_t* b) {
  rational_add(sum, a, b);
}

void lp_rational_add_integer(lp_rational_t* sum, const lp_rational_t* a, const lp_integer_t* b) {
  rational_add_integer(sum, a, b);
}

void lp_rational_sub(lp_rational_t* sub, const lp_rational_t* a, const lp_rational_t* b) {
  rational_sub(sub, a, b);
}

void lp_rational_neg(lp_rational_t* neg, const lp_rational_t* a) {
  rational_neg(neg, a);
}

void lp_rational_inv(lp_rational_t* inv, const lp_rational_t* a) {
  rational_inv(inv, a);
}

void lp_rational_mul(lp_rational_t* mul, const lp_rational_t* a, const lp_rational_t* b) {
  rational_mul(mul, a, b);
}

void lp_rational_mul_2exp(lp_rational_t* mul, const lp_rational_t* a, unsigned n) {
  rational_mul_2exp(mul, a, n);
}

void lp_rational_pow(lp_rational_t* pow, const lp_rational_t*a, unsigned n) {
  rational_pow(pow, a, n);
}

void lp_rational_div(lp_rational_t* div, const lp_rational_t* a, const lp_rational_t* b) {
  rational_div(div, a, b);
}

void lp_rational_div_2exp(lp_rational_t* div, const lp_rational_t* a, unsigned n) {
  rational_div_2exp(div, a, n);
}

void lp_rational_get_num(const lp_rational_t* q, lp_integer_t* num) {
  rational_get_num(q, num);
}

void lp_rational_get_den(const lp_rational_t* q, lp_integer_t* den) {
  rational_get_den(q, den);
}
