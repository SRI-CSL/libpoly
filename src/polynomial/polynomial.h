/*
 * polynomial_internal.h
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#pragma once

#include <polynomial.h>

#include "polynomial/coefficient.h"

struct lp_polynomial_struct {
  /** The actual polynomial representation (so we can use it as a coefficient) */
  coefficient_t data;
  /** Is this an external polynomial (needs checks on function entry) */
  char external;
  /** Context of the polynomial */
  const lp_polynomial_context_t* ctx;
};

void polynomial_construct(lp_polynomial_t* A, const lp_polynomial_context_t* ctx);

void polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from);

void polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from);

void polynomial_construct_simple(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const lp_integer_t* c, lp_variable_t x, unsigned n);

void polynomial_destruct(lp_polynomial_t* A);

lp_polynomial_t* polynomial_alloc(void);

lp_polynomial_t* polynomial_new(const lp_polynomial_context_t* ctx);

void polynomial_set_external(lp_polynomial_t* A);

void polynomial_swap(lp_polynomial_t* A1, lp_polynomial_t* A2);

void polynomial_assign(lp_polynomial_t* A, const lp_polynomial_t* from);

const lp_polynomial_context_t* polynomial_context(const lp_polynomial_t* A);

lp_variable_t polynomial_top_variable(const lp_polynomial_t* A);

size_t polynomial_degree(const lp_polynomial_t* A);

void polynomial_get_coefficient(lp_polynomial_t* C_p, const lp_polynomial_t* A, size_t k);

void polynomial_reductum(lp_polynomial_t* R, const lp_polynomial_t* A);

void polynomial_reductum_m(lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m);

int polynomial_is_constant(const lp_polynomial_t* A);

int polynomial_is_zero(const lp_polynomial_t* A);

int polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m);

int polynomial_cmp(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int polynomial_cmp_type(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int polynomial_divides(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int polynomial_print(const lp_polynomial_t* A, FILE* out);

char* polynomial_to_string(const lp_polynomial_t* A);

int polynomial_same_var_univariate(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_add(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_sub(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_neg(lp_polynomial_t* N, const lp_polynomial_t* A);

void polynomial_mul(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_shl(lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n);

void polynomial_pow(lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n);

void polynomial_add_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_sub_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_div(lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_rem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_divrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_derivative(lp_polynomial_t* A_d, const lp_polynomial_t* A);

void polynomial_gcd(lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_lcm(lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void polynomial_resultant(lp_polynomial_t* res, const lp_polynomial_t* A, const lp_polynomial_t* B);

void polynomial_psc(lp_polynomial_t** psc, const lp_polynomial_t* A, const lp_polynomial_t* B);

void polynomial_reduce(const lp_polynomial_t* A, const lp_polynomial_t* B, lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R);

void polynomial_factor_square_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size);

