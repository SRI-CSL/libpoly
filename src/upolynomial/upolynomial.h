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
typedef struct upolynomial_struct {

  /** The ring of coefficients */
  int_ring K;
  /** The number of monomials */
  size_t size;
  /** The monomials */
  umonomial_t monomials[0];

} upolynomial_t;


typedef struct upolynomial_factors_struct {

  /** Constant factor */
  integer_t constant;
  /** Number of actual factors */
  size_t size;
  /** Size of the factors array */
  size_t capacity;
  /** The irreducible factors */
  upolynomial_t** factors;
  /** The multiplicity of individual factors */
  size_t* multiplicities;

} upolynomial_factors_t;

size_t upolynomial_degree(const upolynomial_t* p);

upolynomial_t* upolynomial_construct_empty(int_ring K, size_t size);

upolynomial_t* upolynomial_construct(int_ring K, size_t degree, const integer_t* coefficients);

void upolynomial_destruct(upolynomial_t* p);

upolynomial_t* upolynomial_construct_power(int_ring K, size_t degree, long c);

upolynomial_t* upolynomial_construct_from_int(int_ring K, size_t degree, const int* coefficients);

upolynomial_t* upolynomial_construct_from_long(int_ring K, size_t degree, const long* coefficients);

upolynomial_t* upolynomial_construct_copy(const upolynomial_t* p);

upolynomial_t* upolynomial_construct_copy_K(int_ring K, const upolynomial_t* p);

void upolynomial_unpack(const upolynomial_t* p, integer_t* out);

int_ring upolynomial_ring(const upolynomial_t* p);

void upolynomial_set_ring(upolynomial_t* p, int_ring K);

const integer_t* upolynomial_lead_coeff(const upolynomial_t* p);

int upolynomial_print(const upolynomial_t* p, FILE* out);

char* upolynomial_to_string(const upolynomial_t* p);

int upolynomial_cmp(const upolynomial_t* p, const upolynomial_t* q);

int upolynomial_is_zero(const upolynomial_t* p);

int upolynomial_is_one(const upolynomial_t* p);

int upolynomial_is_monic(const upolynomial_t* p);

upolynomial_t* upolynomial_add(const upolynomial_t* p, const upolynomial_t* q);

upolynomial_t* upolynomial_sub(const upolynomial_t* p, const upolynomial_t* q);

upolynomial_t* upolynomial_multiply_simple(const umonomial_t* m, const upolynomial_t* q);

upolynomial_t* upolynomial_multiply(const upolynomial_t* p, const upolynomial_t* q);

upolynomial_t* upolynomial_multiply_c(const upolynomial_t* p, const integer_t* c);

upolynomial_t* upolynomial_pow(const upolynomial_t* p, long pow);

upolynomial_t* upolynomial_derivative(const upolynomial_t* p);

upolynomial_t* upolynomial_div_degrees(const upolynomial_t* p, size_t a);

upolynomial_t* upolynomial_div_exact(const upolynomial_t* p, const upolynomial_t* q);

upolynomial_t* upolynomial_div_exact_c(const upolynomial_t* p, const integer_t* c);

upolynomial_t* upolynomial_rem_exact(const upolynomial_t* p, const upolynomial_t* q);

void upolynomial_div_rem_exact(const upolynomial_t* p, const upolynomial_t* q, upolynomial_t** div, upolynomial_t** rem);

void upolynomial_div_pseudo(upolynomial_t** div, upolynomial_t** rem, const upolynomial_t* p, const upolynomial_t* q);

int upolynomial_divides(const upolynomial_t* p, const upolynomial_t* q);

void upolynomial_content_Z(const upolynomial_t* p, integer_t* content);

int upolynomial_is_primitive(const upolynomial_t* p);

void upolynomial_make_primitive_Z(upolynomial_t* p);

upolynomial_t* upolynomial_primitive_part_Z(const upolynomial_t* p);

void upolynomial_evaluate_at_integer(const upolynomial_t* p, const integer_t* x, integer_t* value);

void upolynomial_evaluate_at_dyadic_rational(const upolynomial_t* p, const dyadic_rational_t* x, dyadic_rational_t* value);

int upolynomial_sgn_at_integer(const upolynomial_t* p, const integer_t* x);

int upolynomial_sgn_at_rational(const upolynomial_t* p, const rational_t* x);

int upolynomial_sgn_at_dyadic_rational(const upolynomial_t* p, const dyadic_rational_t* x);

upolynomial_t* upolynomial_gcd(const upolynomial_t* p, const upolynomial_t* q);

upolynomial_t* upolynomial_extended_gcd(const upolynomial_t* p, const upolynomial_t* q, upolynomial_t** u, upolynomial_t** v);

void upolynomial_solve_bezout(const upolynomial_t* p, const upolynomial_t* q, const upolynomial_t* r, upolynomial_t** u, upolynomial_t** v);

upolynomial_factors_t* upolynomial_factor(const upolynomial_t* p);

int upolynomial_roots_count(const upolynomial_t* p, const interval_t* ab);

void upolynomial_roots_isolate(const upolynomial_t* p, algebraic_number_t* roots, size_t* roots_size);

void upolynomial_roots_sturm_sequence(const upolynomial_t* f, upolynomial_t*** S, size_t* size);

upolynomial_factors_t* factors_construct(void);

void factors_swap(upolynomial_factors_t* f1, upolynomial_factors_t* f2);

void factors_clear(upolynomial_factors_t* f);

void factors_destruct(upolynomial_factors_t* f, int destruct_factors);

size_t factors_size(const upolynomial_factors_t* f);

upolynomial_t* factors_get_factor(upolynomial_factors_t* f, size_t i, size_t* d);

const integer_t* factors_get_constant(const upolynomial_factors_t* f);

void factors_add(upolynomial_factors_t* f, upolynomial_t* p, size_t d);

int factors_print(const upolynomial_factors_t* f, FILE* out);

int_ring factors_ring(const upolynomial_factors_t* f);

void factors_set_ring(upolynomial_factors_t* f, int_ring K);

