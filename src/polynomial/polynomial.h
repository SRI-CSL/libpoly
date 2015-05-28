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
  /** Hash */
  size_t hash;
  /** Is this an external polynomial (needs checks on function entry) */
  char external;
  /** Context of the polynomial */
  const lp_polynomial_context_t* ctx;
};

void lp_polynomial_construct(lp_polynomial_t* A, const lp_polynomial_context_t* ctx);

void lp_polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from);

void lp_polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from);

void lp_polynomial_construct_simple(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const lp_integer_t* c, lp_variable_t x, unsigned n);

void lp_polynomial_destruct(lp_polynomial_t* A);

lp_polynomial_t* lp_polynomial_alloc(void);

lp_polynomial_t* lp_polynomial_new(const lp_polynomial_context_t* ctx);

void lp_polynomial_set_external(lp_polynomial_t* A);

void lp_polynomial_swap(lp_polynomial_t* A1, lp_polynomial_t* A2);

void lp_polynomial_assign(lp_polynomial_t* A, const lp_polynomial_t* from);

const lp_polynomial_context_t* lp_polynomial_context(const lp_polynomial_t* A);

lp_variable_t lp_polynomial_top_variable(const lp_polynomial_t* A);

size_t lp_polynomial_degree(const lp_polynomial_t* A);

void lp_polynomial_get_coefficient(lp_polynomial_t* C_p, const lp_polynomial_t* A, size_t k);

void lp_polynomial_reductum(lp_polynomial_t* R, const lp_polynomial_t* A);

void lp_polynomial_reductum_m(lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m);

int lp_polynomial_is_constant(const lp_polynomial_t* A);

int lp_polynomial_is_zero(const lp_polynomial_t* A);

int lp_polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m);

int lp_polynomial_cmp(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int lp_polynomial_cmp_type(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int lp_polynomial_divides(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

int lp_polynomial_print(const lp_polynomial_t* A, FILE* out);

char* lp_polynomial_to_string(const lp_polynomial_t* A);

int polynomial_same_var_univariate(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_add(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_sub(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_neg(lp_polynomial_t* N, const lp_polynomial_t* A);

void lp_polynomial_mul(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_shl(lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n);

void lp_polynomial_pow(lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n);

void lp_polynomial_add_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_sub_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_div(lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_rem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_divrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_derivative(lp_polynomial_t* A_d, const lp_polynomial_t* A);

void lp_polynomial_gcd(lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_lcm(lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

void lp_polynomial_resultant(lp_polynomial_t* res, const lp_polynomial_t* A, const lp_polynomial_t* B);

void lp_polynomial_psc(lp_polynomial_t** psc, const lp_polynomial_t* A, const lp_polynomial_t* B);

void lp_polynomial_reduce(const lp_polynomial_t* A, const lp_polynomial_t* B, lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R);

void lp_polynomial_factor_square_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size);

