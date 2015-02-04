/*
 * dyadic_dyadic_rational.c
 *
 *  Created on: Jan 22, 2014
 *      Author: dejan
 */

#include "number/dyadic_rational.h"

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
