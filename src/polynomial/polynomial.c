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

void lp_polynomial_external_clean(const lp_polynomial_t* A_const) {
  if (A_const->external && !coefficient_in_order(A_const->ctx, &A_const->data)) {
    lp_polynomial_t* A = (lp_polynomial_t*) A_const;
    coefficient_order(A->ctx, &A->data);
  }
}

void lp_polynomial_set_context(lp_polynomial_t* A, const lp_polynomial_context_t* ctx) {
  if (A->ctx != ctx) {
    if (A->ctx && A->external) {
      lp_polynomial_context_detach((lp_polynomial_context_t*)A->ctx);
    }
    A->ctx = ctx;
    if (A->ctx && A->external) {
      lp_polynomial_context_attach((lp_polynomial_context_t*)A->ctx);
    }
  }
}

void lp_polynomial_construct(lp_polynomial_t* A, const lp_polynomial_context_t* ctx) {
  A->ctx = 0;
  A->external = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct(ctx, &A->data);
}

void lp_polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from) {
  A->ctx = 0;
  A->external = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct_copy(A->ctx, &A->data, from);
}

void lp_polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from) {
  A->ctx = 0;
  A->external = 0;
  lp_polynomial_set_context(A, from->ctx);
  coefficient_construct_copy(A->ctx, &A->data, &from->data);
}

/** Construct a simple polynomial c*x^n */
void lp_polynomial_construct_simple(
    lp_polynomial_t* A, const lp_polynomial_context_t* ctx,
    const lp_integer_t* c, lp_variable_t x, unsigned n)
{
  A->ctx = 0;
  A->external = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct_simple(ctx, &A->data, c, x, n);
}

void lp_polynomial_destruct(lp_polynomial_t* A) {
  coefficient_destruct(&A->data);
  if (A->external) {
    lp_polynomial_context_detach((lp_polynomial_context_t*)A->ctx);
  }
}

lp_polynomial_t* lp_polynomial_alloc(void) {
  lp_polynomial_t* new = malloc(sizeof(lp_polynomial_t));
  return new;
}

lp_polynomial_t* lp_polynomial_new(const lp_polynomial_context_t* ctx) {
  lp_polynomial_t* new = lp_polynomial_alloc();
  lp_polynomial_construct(new, ctx);
  return new;
}

void lp_polynomial_set_external(lp_polynomial_t* A) {
  if (!A->external) {
    A->external = 1;
    lp_polynomial_context_attach((lp_polynomial_context_t*) A->ctx);
  }
}

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

void lp_polynomial_swap(lp_polynomial_t* A1, lp_polynomial_t* A2) {
  // Swap everything, but keep the external flags
  lp_polynomial_t tmp = *A1; *A1 = *A2; *A2 = tmp;
  SWAP(unsigned, A1->external, A2->external);
}

void lp_polynomial_assign(lp_polynomial_t* A, const lp_polynomial_t* from) {
  if (A != from) {
    lp_polynomial_set_context(A, from->ctx);
    coefficient_assign(A->ctx, &A->data, &from->data);
  }
}

const lp_polynomial_context_t* lp_polynomial_context(const lp_polynomial_t* A) {
  return A->ctx;
}

lp_variable_t lp_polynomial_top_variable(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_top_variable(&A->data);
}

size_t lp_polynomial_degree(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_degree(&A->data);
}

void lp_polynomial_get_coefficient(lp_polynomial_t* C_p, const lp_polynomial_t* A, size_t k) {
  lp_polynomial_external_clean(A);

  if (k > lp_polynomial_degree(A)) {
    lp_polynomial_t result;
    lp_polynomial_construct(&result, A->ctx);
    lp_polynomial_swap(C_p, &result);
    lp_polynomial_destruct(&result);
  } else {
    const coefficient_t* C = coefficient_get_coefficient(&A->data, k);
    lp_polynomial_t result;
    lp_polynomial_construct_from_coefficient(&result, A->ctx, C);
    lp_polynomial_swap(C_p, &result);
    lp_polynomial_destruct(&result);
  }
}

void lp_polynomial_reductum(lp_polynomial_t* R, const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  lp_polynomial_set_context(R, A->ctx);
  coefficient_reductum(A->ctx, &R->data, &A->data);
}

void lp_polynomial_reductum_m(lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);
  lp_polynomial_set_context(R, A->ctx);
  coefficient_reductum_m(A->ctx, &R->data, &A->data, m);
}

int lp_polynomial_is_constant(const lp_polynomial_t* A) {
  return coefficient_is_constant(&A->data);
}

int lp_polynomial_is_zero(const lp_polynomial_t* A) {
  return coefficient_is_zero(A->ctx, &A->data);
}

int lp_polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m) {
  return coefficient_sgn(A->ctx, &A->data, m);
}

int lp_polynomial_cmp(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_cmp("); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
  }

  if (!lp_polynomial_context_equal(A1->ctx, A2->ctx)) {
    // random order for different contexts
    return A1 - A2;
  }

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);
  int cmp = coefficient_cmp(A1->ctx, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_cmp("); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(") => %d\n", cmp);
  }

  return cmp;
}

int lp_polynomial_cmp_type(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;
  assert(lp_polynomial_context_equal(A1->ctx, ctx));
  assert(lp_polynomial_context_equal(A2->ctx, ctx));
  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);
  return coefficient_cmp_type(ctx, &A1->data, &A2->data);
}

int lp_polynomial_divides(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  if (!lp_polynomial_context_equal(A1->ctx, A2->ctx)) {
    return 0;
  }
  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);
  return coefficient_divides(A1->ctx, &A1->data, &A2->data);
}

int lp_polynomial_print(const lp_polynomial_t* A, FILE* out) {
  lp_polynomial_external_clean(A);
  return coefficient_print(A->ctx, &A->data, out);
}

char* lp_polynomial_to_string(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_to_string(A->ctx, &A->data);
}

void lp_polynomial_add(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_add("); lp_polynomial_print(S, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(S, A1->ctx);

  coefficient_add(S->ctx, &S->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_add() => "); lp_polynomial_print(S, trace_out); tracef("\n");
  }
}

void lp_polynomial_sub(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sub("); lp_polynomial_print(S, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(S, A1->ctx);

  coefficient_sub(S->ctx, &S->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sub() => "); lp_polynomial_print(S, trace_out); tracef("\n");
  }
}

void lp_polynomial_neg(lp_polynomial_t* N, const lp_polynomial_t* A) {

  lp_polynomial_external_clean(A);

  lp_polynomial_set_context(N, N->ctx);

  coefficient_neg(N->ctx, &N->data, &A->data);
}

void lp_polynomial_mul(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mul("); lp_polynomial_print(P, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(P, A1->ctx);

  coefficient_mul(P->ctx, &P->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mul() => "); lp_polynomial_print(P, trace_out); tracef("\n");
  }
}

void lp_polynomial_shl(lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n) {

  lp_polynomial_external_clean(A);

  lp_polynomial_set_context(S, A->ctx);

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  coefficient_shl(S->ctx, &S->data, &A->data, VAR(&A->data), n);
}

void lp_polynomial_pow(lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_pow("); lp_polynomial_print(P, trace_out); tracef(", "); lp_polynomial_print(A, trace_out); tracef(")\n");
    lp_variable_order_print(
        A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);

  lp_polynomial_set_context(P, A->ctx);

  coefficient_pow(P->ctx, &P->data, &A->data, n);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_pow() => "); lp_polynomial_print(P, trace_out); tracef("\n");
  }
}

void lp_polynomial_add_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;

  assert(lp_polynomial_context_equal(S->ctx, ctx));
  assert(lp_polynomial_context_equal(A1->ctx, ctx));
  assert(lp_polynomial_context_equal(A2->ctx, ctx));

  lp_polynomial_external_clean(S);
  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  coefficient_add_mul(ctx, &S->data, &A1->data, &A2->data);
}

void lp_polynomial_sub_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  const lp_polynomial_context_t* ctx = A1->ctx;

  assert(lp_polynomial_context_equal(S->ctx, ctx));
  assert(lp_polynomial_context_equal(A1->ctx, ctx));
  assert(lp_polynomial_context_equal(A2->ctx, ctx));

  lp_polynomial_external_clean(S);
  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  coefficient_sub_mul(ctx, &S->data, &A1->data, &A2->data);
}

void lp_polynomial_div(lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_div("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(D, A1->ctx);

  coefficient_div(D->ctx, &D->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_div() => "); lp_polynomial_print(D, trace_out); tracef("\n");
  }
}

void lp_polynomial_rem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem("); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(R, A1->ctx);

  coefficient_rem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem() => "); lp_polynomial_print(R, trace_out); tracef("\n");
  }
}

void lp_polynomial_prem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_prem("); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(R, A1->ctx);

  coefficient_prem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_prem() => "); lp_polynomial_print(R, trace_out); tracef("\n");
  }
}

void lp_polynomial_sprem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sprem("); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(R, A1->ctx);

  coefficient_sprem(R->ctx, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_sprem() => "); lp_polynomial_print(R, trace_out); tracef("\n");
  }
}

void lp_polynomial_divrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_divrem("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(D, A1->ctx);
  lp_polynomial_set_context(R, A1->ctx);

  coefficient_divrem(D->ctx, &R->data, &R->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_rem() => ("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(")\n");
  }
}

void lp_polynomial_derivative(lp_polynomial_t* A_d, const lp_polynomial_t* A) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative("); lp_polynomial_print(A_d, trace_out); tracef(", "); lp_polynomial_print(A, trace_out); tracef(")\n");
    lp_variable_order_print(
        A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);

  lp_polynomial_set_context(A_d, A->ctx);

  coefficient_derivative(A_d->ctx, &A_d->data, &A->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative() => "); lp_polynomial_print(A_d, trace_out); tracef("\n");
  }
}

void lp_polynomial_gcd(lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_gcd("); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(gcd, A1->ctx);

  coefficient_gcd(gcd->ctx, &gcd->data, &A1->data, &A2->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_gcd() => "); lp_polynomial_print(gcd, trace_out); tracef("\n");
  }
}

void lp_polynomial_lcm(lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {
  assert(lp_polynomial_context_equal(A1->ctx, A2->ctx));

  lp_polynomial_external_clean(A1);
  lp_polynomial_external_clean(A2);

  lp_polynomial_set_context(lcm, A1->ctx);

  coefficient_lcm(lcm->ctx, &lcm->data, &A1->data, &A2->data);
}

void lp_polynomial_reduce(
    const lp_polynomial_t* A, const lp_polynomial_t* B,
    lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R)
{
  const lp_polynomial_context_t* ctx = A->ctx;

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_reduce("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(")\n");
    lp_variable_order_print(
        A->ctx->var_order, A->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  assert(lp_polynomial_context_equal(B->ctx, ctx));

  lp_polynomial_external_clean(A);
  lp_polynomial_external_clean(B);

  lp_polynomial_set_context(P, ctx);
  lp_polynomial_set_context(Q, ctx);
  lp_polynomial_set_context(R, ctx);

  coefficient_reduce(ctx, &A->data, &B->data, &P->data, &Q->data, &R->data, 1);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_derivative() =>\n");
    tracef("\t P = "); lp_polynomial_print(P, trace_out); tracef("\n");
    tracef("\t Q = "); lp_polynomial_print(Q, trace_out); tracef("\n");
    tracef("\t R = "); lp_polynomial_print(R, trace_out); tracef("\n");
  }
}

void lp_polynomial_psc(lp_polynomial_t** psc, const lp_polynomial_t* A, const lp_polynomial_t* B) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_psc("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(")\n");
  }

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  size_t A_deg = lp_polynomial_degree(A);
  size_t B_deg = lp_polynomial_degree(B);

  if (A_deg < B_deg) {
    lp_polynomial_psc(psc, B, A);
    return;
  }

  const lp_polynomial_context_t* ctx = A->ctx;
  assert(lp_polynomial_context_equal(B->ctx, ctx));

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);
  lp_polynomial_external_clean(B);

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
    lp_polynomial_construct_from_coefficient(&tmp, ctx, psc_coeff + i);
    lp_polynomial_swap(&tmp, psc[i]);
    lp_polynomial_destruct(&tmp);
    coefficient_destruct(&psc_coeff[i]);
  }

  free(psc_coeff);

  if (trace_is_enabled("polynomial")) {
    for (i = 0; i < size; ++ i) {
      tracef("PSC[%zu] = ", i); lp_polynomial_print(psc[i], trace_out); tracef("\n");
    }
  }
}

void lp_polynomial_resultant(lp_polynomial_t* res, const lp_polynomial_t* A, const lp_polynomial_t* B) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_resultant("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(")\n");
  }

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  const lp_polynomial_context_t* ctx = A->ctx;
  assert(lp_polynomial_context_equal(B->ctx, ctx));

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);
  lp_polynomial_external_clean(B);

  // Compute
  coefficient_resultant(ctx, &res->data, &A->data, &B->data);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_resultant("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(") => "); lp_polynomial_print(res, trace_out); tracef("\n");
  }
}

void lp_polynomial_factor_square_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_factor_square_free("); lp_polynomial_print(A, trace_out); tracef(")\n");
  }

  assert(*factors == 0);
  assert(*multiplicities == 0);
  assert(*size == 0);

  const lp_polynomial_context_t* ctx = A->ctx;

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);
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
    lp_polynomial_construct_from_coefficient((*factors)[i], A->ctx, coeff_factors.factors + i);
    (*multiplicities)[i] = coeff_factors.multiplicities[i];
  }

  coefficient_factors_destruct(&coeff_factors);
}
