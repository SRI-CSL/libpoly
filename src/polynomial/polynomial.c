/*
 * polynomial_internal.c
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#include "polynomial/polynomial.h"

#include "polynomial/gcd.h"
#include "polynomial/factorization.h"
#include "polynomial/output.h"

#include "utils/debug_trace.h"

#include <assert.h>
#include <malloc.h>

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

void polynomial_external_clean(const lp_polynomial_t* A_const) {
  if (A_const->external && !coefficient_in_order(A_const->ctx, &A_const->data)) {
    lp_polynomial_t* A = (lp_polynomial_t*) A_const;
    coefficient_order(A->ctx, &A->data);
  }
}

void polynomial_set_context(lp_polynomial_t* A, const lp_polynomial_context_t* ctx) {
  if (A->ctx != ctx) {
    if (A->ctx && A->external) {
      polynomial_context_ops.detach((lp_polynomial_context_t*)A->ctx);
    }
    A->ctx = ctx;
    if (A->ctx && A->external) {
      polynomial_context_ops.attach((lp_polynomial_context_t*)A->ctx);
    }
  }
}

void polynomial_construct(lp_polynomial_t* A, const lp_polynomial_context_t* ctx) {
  A->ctx = 0;
  A->external = 0;
  polynomial_set_context(A, ctx);
  coefficient_construct(ctx, &A->data);
}

void polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from) {
  A->ctx = 0;
  A->external = 0;
  polynomial_set_context(A, ctx);
  coefficient_construct_copy(A->ctx, &A->data, from);
}

void polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from) {
  A->ctx = 0;
  A->external = 0;
  polynomial_set_context(A, from->ctx);
  coefficient_construct_copy(A->ctx, &A->data, &from->data);
}

/** Construct a simple polynomial c*x^n */
void polynomial_construct_simple(
    lp_polynomial_t* A, const lp_polynomial_context_t* ctx,
    const lp_integer_t* c, lp_variable_t x, unsigned n)
{
  A->ctx = 0;
  A->external = 0;
  polynomial_set_context(A, ctx);
  coefficient_construct_simple(ctx, &A->data, c, x, n);
}

void polynomial_destruct(lp_polynomial_t* A) {
  coefficient_destruct(&A->data);
  if (A->external) {
    polynomial_context_ops.detach((lp_polynomial_context_t*)A->ctx);
  }
}

lp_polynomial_t* polynomial_alloc(void) {
  lp_polynomial_t* new = malloc(sizeof(lp_polynomial_t));
  return new;
}

lp_polynomial_t* polynomial_new(const lp_polynomial_context_t* ctx) {
  lp_polynomial_t* new = polynomial_alloc();
  polynomial_construct(new, ctx);
  return new;
}

void polynomial_set_external(lp_polynomial_t* A) {
  if (!A->external) {
    A->external = 1;
    polynomial_context_ops.attach((lp_polynomial_context_t*) A->ctx);
  }
}

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

void polynomial_swap(lp_polynomial_t* A1, lp_polynomial_t* A2) {
  // Swap everything, but keep the external flags
  lp_polynomial_t tmp = *A1; *A1 = *A2; *A2 = tmp;
  SWAP(unsigned, A1->external, A2->external);
}

void polynomial_assign(lp_polynomial_t* A, const lp_polynomial_t* from) {
  if (A != from) {
    polynomial_set_context(A, from->ctx);
    coefficient_assign(A->ctx, &A->data, &from->data);
  }
}

const lp_polynomial_context_t* polynomial_context(const lp_polynomial_t* A) {
  return A->ctx;
}

lp_variable_t polynomial_top_variable(const lp_polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_top_variable(&A->data);
}

size_t polynomial_degree(const lp_polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_degree(&A->data);
}

void polynomial_get_coefficient(lp_polynomial_t* C_p, const lp_polynomial_t* A, size_t k) {
  polynomial_external_clean(A);

  if (k > polynomial_degree(A)) {
    lp_polynomial_t result;
    polynomial_construct(&result, A->ctx);
    polynomial_swap(C_p, &result);
    polynomial_destruct(&result);
  } else {
    const coefficient_t* C = coefficient_get_coefficient(&A->data, k);
    lp_polynomial_t result;
    polynomial_construct_from_coefficient(&result, A->ctx, C);
    polynomial_swap(C_p, &result);
    polynomial_destruct(&result);
  }
}

void polynomial_reductum(lp_polynomial_t* R, const lp_polynomial_t* A) {
  polynomial_external_clean(A);
  polynomial_set_context(R, A->ctx);
  coefficient_reductum(A->ctx, &R->data, &A->data);
}

void polynomial_reductum_m(lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m) {
  polynomial_external_clean(A);
  polynomial_set_context(R, A->ctx);
  coefficient_reductum_m(A->ctx, &R->data, &A->data, m);
}

int polynomial_is_constant(const lp_polynomial_t* A) {
  return coefficient_is_constant(&A->data);
}

int polynomial_is_zero(const lp_polynomial_t* A) {
  return coefficient_is_zero(A->ctx, &A->data);
}

int polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m) {
  return coefficient_sgn(A->ctx, &A->data, m);
}

int polynomial_cmp(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_cmp("); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
  }

  if (!polynomial_context_ops.equal(A1->ctx, A2->ctx)) {
    // random order for different contexts
    return A1 - A2;
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  int cmp = coefficient_cmp(A1->ctx, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_cmp("); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(") => %d\n", cmp);
  }

  return cmp;
}

int polynomial_cmp_type(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  return coefficient_cmp_type(ctx, &A1->data, &A2->data);
}

int polynomial_divides(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  if (!polynomial_context_ops.equal(A1->ctx, A2->ctx)) {
    return 0;
  }
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  return coefficient_divides(A1->ctx, &A1->data, &A2->data);
}

int polynomial_print(const lp_polynomial_t* A, FILE* out) {
  polynomial_external_clean(A);
  return coefficient_print(A->ctx, &A->data, out);
}

char* polynomial_to_string(const lp_polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_to_string(A->ctx, &A->data);
}

void polynomial_add(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_add("); polynomial_print(S, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(S, A1->ctx);

  coefficient_add(S->ctx, &S->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_add() => "); polynomial_print(S, trace_out); tracef("\n");
  }
}

void polynomial_sub(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sub("); polynomial_print(S, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(S, A1->ctx);

  coefficient_sub(S->ctx, &S->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sub() => "); polynomial_print(S, trace_out); tracef("\n");
  }
}

void polynomial_neg(lp_polynomial_t* N, const lp_polynomial_t* A) {

  polynomial_external_clean(A);

  polynomial_set_context(N, N->ctx);

  coefficient_neg(N->ctx, &N->data, &A->data);
}

void polynomial_mul(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mul("); polynomial_print(P, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(P, A1->ctx);

  coefficient_mul(P->ctx, &P->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mul() => "); polynomial_print(P, trace_out); tracef("\n");
  }
}

void polynomial_shl(lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n) {

  polynomial_external_clean(A);

  polynomial_set_context(S, A->ctx);

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  coefficient_shl(S->ctx, &S->data, &A->data, VAR(&A->data), n);
}

void polynomial_pow(lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_pow("); polynomial_print(P, trace_out); tracef(", "); polynomial_print(A, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);

  polynomial_set_context(P, A->ctx);

  coefficient_pow(P->ctx, &P->data, &A->data, n);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_pow() => "); polynomial_print(P, trace_out); tracef("\n");
  }
}

void polynomial_add_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;

  assert(polynomial_context_ops.equal(S->ctx, ctx));
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));

  polynomial_external_clean(S);
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  coefficient_add_mul(ctx, &S->data, &A1->data, &A2->data);
}

void polynomial_sub_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;

  assert(polynomial_context_ops.equal(S->ctx, ctx));
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));

  polynomial_external_clean(S);
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  coefficient_sub_mul(ctx, &S->data, &A1->data, &A2->data);
}

void polynomial_div(lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_div("); polynomial_print(D, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(D, A1->ctx);

  coefficient_div(D->ctx, &D->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_div() => "); polynomial_print(D, trace_out); tracef("\n");
  }
}

void polynomial_rem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem("); polynomial_print(R, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(R, A1->ctx);

  coefficient_rem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem() => "); polynomial_print(R, trace_out); tracef("\n");
  }
}

void polynomial_prem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_prem("); polynomial_print(R, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(R, A1->ctx);

  coefficient_prem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_prem() => "); polynomial_print(R, trace_out); tracef("\n");
  }
}

void polynomial_sprem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sprem("); polynomial_print(R, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(R, A1->ctx);

  coefficient_sprem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sprem() => "); polynomial_print(R, trace_out); tracef("\n");
  }
}

void polynomial_divrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_divrem("); polynomial_print(D, trace_out); tracef(", "); polynomial_print(R, trace_out); tracef(", "); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(D, A1->ctx);
  polynomial_set_context(R, A1->ctx);

  coefficient_divrem(D->ctx, &R->data, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem() => ("); polynomial_print(D, trace_out); tracef(", "); polynomial_print(R, trace_out); tracef(")\n");
  }
}

void polynomial_derivative(lp_polynomial_t* A_d, const lp_polynomial_t* A) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative("); polynomial_print(A_d, trace_out); tracef(", "); polynomial_print(A, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);

  polynomial_set_context(A_d, A->ctx);

  coefficient_derivative(A_d->ctx, &A_d->data, &A->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative() => "); polynomial_print(A_d, trace_out); tracef("\n");
  }
}

void polynomial_gcd(lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_gcd("); polynomial_print(A1, trace_out); tracef(", "); polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(gcd, A1->ctx);

  coefficient_gcd(gcd->ctx, &gcd->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_gcd() => "); polynomial_print(gcd, trace_out); tracef("\n");
  }
}

void polynomial_lcm(lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(lcm, A1->ctx);

  coefficient_lcm(lcm->ctx, &lcm->data, &A1->data, &A2->data);
}

void polynomial_reduce(
    const lp_polynomial_t* A, const lp_polynomial_t* B,
    lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R)
{
  const lp_polynomial_context_t* ctx = A->ctx;

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_reduce("); polynomial_print(A, trace_out); tracef(", "); polynomial_print(B, trace_out); tracef(")\n");
    lp_variable_order_simple_ops.print(
        (lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(polynomial_context_ops.equal(B->ctx, ctx));

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  polynomial_set_context(P, ctx);
  polynomial_set_context(Q, ctx);
  polynomial_set_context(R, ctx);

  coefficient_reduce(ctx, &A->data, &B->data, &P->data, &Q->data, &R->data, 1);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative() =>\n");
    tracef("\t P = "); polynomial_print(P, trace_out); tracef("\n");
    tracef("\t Q = "); polynomial_print(Q, trace_out); tracef("\n");
    tracef("\t R = "); polynomial_print(R, trace_out); tracef("\n");
  }
}

void polynomial_psc(lp_polynomial_t** psc, const lp_polynomial_t* A, const lp_polynomial_t* B) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_psc("); polynomial_print(A, trace_out); tracef(", "); polynomial_print(B, trace_out); tracef(")\n");
  }

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  size_t A_deg = polynomial_degree(A);
  size_t B_deg = polynomial_degree(B);

  if (A_deg < B_deg) {
    polynomial_psc(psc, B, A);
    return;
  }

  const lp_polynomial_context_t* ctx = A->ctx;
  assert(polynomial_context_ops.equal(B->ctx, ctx));

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_simple_ops.print((lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  // Allocate the space for the result
  size_t size = B_deg + 1;
  coefficient_t* psc_coeff = malloc(sizeof(coefficient_t)*size);
  size_t i;
  for (i = 0; i < size; ++ i) {
    coefficient_construct(ctx, psc_coeff + i);
  }

  // Compute
  coefficient_psc(ctx, psc_coeff, &A->data, &B->data);

  // Construct the output (one less, we ignore the final 1)
  for (i = 0; i < size; ++ i) {
    lp_polynomial_t tmp;
    polynomial_construct_from_coefficient(&tmp, ctx, psc_coeff + i);
    polynomial_swap(&tmp, psc[i]);
    polynomial_destruct(&tmp);
    coefficient_destruct(&psc_coeff[i]);
  }

  free(psc_coeff);

  if (trace_is_enabled("polynomial")) {
    for (i = 0; i < size; ++ i) {
      tracef("PSC[%zu] = ", i); polynomial_print(psc[i], trace_out); tracef("\n");
    }
  }
}

void polynomial_resultant(lp_polynomial_t* res, const lp_polynomial_t* A, const lp_polynomial_t* B) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_resultant("); polynomial_print(A, trace_out); tracef(", "); polynomial_print(B, trace_out); tracef(")\n");
  }

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  const lp_polynomial_context_t* ctx = A->ctx;
  assert(polynomial_context_ops.equal(B->ctx, ctx));

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_simple_ops.print((lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  // Compute
  coefficient_resultant(ctx, &res->data, &A->data, &B->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_resultant("); polynomial_print(A, trace_out); tracef(", "); polynomial_print(B, trace_out); tracef(") => "); polynomial_print(res, trace_out); tracef("\n");
  }
}

void polynomial_factor_square_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_factor_square_free("); polynomial_print(A, trace_out); tracef(")\n");
  }

  assert(*factors == 0);
  assert(*multiplicities == 0);
  assert(*size == 0);

  const lp_polynomial_context_t* ctx = A->ctx;

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_simple_ops.print((lp_variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  coefficient_factors_t coeff_factors;
  coefficient_factors_construct(&coeff_factors);

  coefficient_factor_square_free(ctx, &A->data, &coeff_factors);

  if (coeff_factors.size) {
    *size = coeff_factors.size;
    *factors = malloc(sizeof(lp_polynomial_t*) * (*size));
    *multiplicities = malloc(sizeof(size_t) * (*size));
  } else {
    *size = 0;
    *factors = 0;
    *multiplicities = 0;
  }

  size_t i;
  for (i = 0; i < *size; ++ i) {
    (*factors)[i] = malloc(sizeof(lp_polynomial_t));
    polynomial_construct_from_coefficient((*factors)[i], A->ctx, coeff_factors.factors + i);
    (*multiplicities)[i] = coeff_factors.multiplicities[i];
  }

  coefficient_factors_destruct(&coeff_factors);
}

const polynomial_ops_t polynomial_ops = {
  polynomial_construct,
  polynomial_construct_simple,
  polynomial_construct_copy,
  polynomial_destruct,
  polynomial_alloc,
  polynomial_new,
  polynomial_set_external,
  polynomial_swap,
  polynomial_assign,
  polynomial_context,
  polynomial_degree,
  polynomial_top_variable,
  polynomial_get_coefficient,
  polynomial_reductum,
  polynomial_reductum_m,
  polynomial_is_constant,
  polynomial_is_zero,
  polynomial_sgn,
  polynomial_cmp,
  polynomial_cmp_type,
  polynomial_divides,
  polynomial_print,
  polynomial_to_string,
  polynomial_add,
  polynomial_sub,
  polynomial_neg,
  polynomial_mul,
  polynomial_shl,
  polynomial_pow,
  polynomial_add_mul,
  polynomial_sub_mul,
  polynomial_reduce,
  polynomial_div,
  polynomial_rem,
  polynomial_prem,
  polynomial_sprem,
  polynomial_divrem,
  polynomial_derivative,
  polynomial_gcd,
  polynomial_lcm,
  polynomial_resultant,
  polynomial_psc,
  polynomial_factor_square_free
};

