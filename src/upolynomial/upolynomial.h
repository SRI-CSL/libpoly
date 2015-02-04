/*
 * univariate_polynomial_internal.h
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#pragma once

#include <upolynomial.h>

#include "umonomial.h"

/**
 * A polynomial is the ring, number of monomials and the monomials.
 */
struct lp_upolynomial_struct {
  /** The ring of coefficients */
  lp_int_ring K;
  /** The number of monomials */
  size_t size;
  /** The monomials */
  umonomial_t monomials[];
};

size_t upolynomial_degree(const lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_construct_empty(lp_int_ring K, size_t size);

lp_upolynomial_t* upolynomial_construct(lp_int_ring K, size_t degree, const lp_integer_t* coefficients);

void upolynomial_delete(lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_construct_power(lp_int_ring K, size_t degree, long c);

lp_upolynomial_t* upolynomial_construct_from_int(lp_int_ring K, size_t degree, const int* coefficients);

lp_upolynomial_t* upolynomial_construct_from_long(lp_int_ring K, size_t degree, const long* coefficients);

lp_upolynomial_t* upolynomial_construct_copy(const lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_construct_copy_K(lp_int_ring K, const lp_upolynomial_t* p);

void upolynomial_unpack(const lp_upolynomial_t* p, lp_integer_t* out);

lp_int_ring upolynomial_ring(const lp_upolynomial_t* p);

void upolynomial_set_ring(lp_upolynomial_t* p, lp_int_ring K);

const lp_integer_t* upolynomial_lead_coeff(const lp_upolynomial_t* p);
\
int upolynomial_cmp(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

int upolynomial_is_zero(const lp_upolynomial_t* p);

int upolynomial_is_one(const lp_upolynomial_t* p);

int upolynomial_is_monic(const lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_add(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

lp_upolynomial_t* upolynomial_sub(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

lp_upolynomial_t* upolynomial_neg(const lp_upolynomial_t* p);

void upolynomial_neg_in_place(lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_mul(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

lp_upolynomial_t* upolynomial_mul_c(const lp_upolynomial_t* p, const lp_integer_t* c);

lp_upolynomial_t* upolynomial_pow(const lp_upolynomial_t* p, long pow);

lp_upolynomial_t* upolynomial_derivative(const lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_div_degrees(const lp_upolynomial_t* p, size_t a);

lp_upolynomial_t* upolynomial_div_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

lp_upolynomial_t* upolynomial_div_exact_c(const lp_upolynomial_t* p, const lp_integer_t* c);

lp_upolynomial_t* upolynomial_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

void upolynomial_div_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q, lp_upolynomial_t** div, lp_upolynomial_t** rem);

void upolynomial_div_pseudo(lp_upolynomial_t** div, lp_upolynomial_t** rem, const lp_upolynomial_t* p, const lp_upolynomial_t* q);

int upolynomial_divides(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

void upolynomial_content_Z(const lp_upolynomial_t* p, lp_integer_t* content);

int upolynomial_is_primitive(const lp_upolynomial_t* p);

void upolynomial_make_primitive_Z(lp_upolynomial_t* p);

lp_upolynomial_t* upolynomial_primitive_part_Z(const lp_upolynomial_t* p);

void upolynomial_evaluate_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x, lp_integer_t* value);

void upolynomial_evaluate_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x, lp_dyadic_rational_t* value);

int upolynomial_sgn_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x);

int upolynomial_sgn_at_rational(const lp_upolynomial_t* p, const lp_rational_t* x);

int upolynomial_sgn_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x);

lp_upolynomial_t* upolynomial_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

lp_upolynomial_t* upolynomial_extended_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q, lp_upolynomial_t** u, lp_upolynomial_t** v);

void upolynomial_solve_bezout(const lp_upolynomial_t* p, const lp_upolynomial_t* q, const lp_upolynomial_t* r, lp_upolynomial_t** u, lp_upolynomial_t** v);

lp_upolynomial_factors_t* upolynomial_factor(const lp_upolynomial_t* p);

int upolynomial_roots_count(const lp_upolynomial_t* p, const lp_interval_t* ab);

void upolynomial_roots_isolate(const lp_upolynomial_t* p, lp_algebraic_number_t* roots, size_t* roots_size);

void upolynomial_roots_sturm_sequence(const lp_upolynomial_t* f, lp_upolynomial_t*** S, size_t* size);

/** Compute f(-x) */
lp_upolynomial_t* upolynomial_subst_x_neg(const lp_upolynomial_t* f);


