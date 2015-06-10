/*
 * polynomial_internal.c
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#include <upolynomial.h>
#include <feasibility_set.h>
#include <variable_db.h>
#include <variable_list.h>

#include "polynomial/polynomial.h"

#include "polynomial/gcd.h"
#include "polynomial/factorization.h"
#include "polynomial/output.h"

#include "number/rational.h"
#include "number/integer.h"

#include "polynomial/feasibility_set.h"

#include "utils/debug_trace.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

static
void check_polynomial_assignment(const lp_polynomial_t* A, const lp_assignment_t* M, lp_variable_t x) {
  assert(A->data.type == COEFFICIENT_NUMERIC || VAR(&A->data) == x);
  lp_variable_list_t vars;
  lp_variable_list_construct(&vars);
  lp_polynomial_get_variables(A, &vars);
  size_t i = 0;
  for (i = 0; i < vars.list_size; ++i) {
    lp_variable_t y = vars.list[i];
    if (x != y && lp_assignment_get_value(M, y)->type == LP_VALUE_NONE) {
      lp_polynomial_print(A, trace_out);
      tracef("\n")
      assert(0);
    }
  }
  lp_variable_list_destruct(&vars);
}

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
  A->hash = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct(ctx, &A->data);
}

void lp_polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from) {
  A->ctx = 0;
  A->external = 0;
  A->hash = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct_copy(A->ctx, &A->data, from);
}

void lp_polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from) {
  A->ctx = 0;
  A->external = 0;
  A->hash = from->hash;
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
  A->hash = 0;
  lp_polynomial_set_context(A, ctx);
  coefficient_construct_simple(ctx, &A->data, c, x, n);
}

void lp_polynomial_destruct(lp_polynomial_t* A) {
  coefficient_destruct(&A->data);
  if (A->external) {
    lp_polynomial_context_detach((lp_polynomial_context_t*)A->ctx);
  }
}

void lp_polynomial_delete(lp_polynomial_t* A) {
  lp_polynomial_destruct(A);
  free(A);
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

lp_polynomial_t* lp_polynomial_new_copy(const lp_polynomial_t* A) {
  lp_polynomial_t* new = lp_polynomial_alloc();
  lp_polynomial_construct_copy(new, A);
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

lp_variable_t lp_polynomial_top_variable(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_top_variable(&A->data);
}

size_t lp_polynomial_degree(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_degree(&A->data);
}

const lp_polynomial_context_t* lp_polynomial_get_context(const lp_polynomial_t* A) {
  return A->ctx;
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

int lp_polynomial_is_univariate(const lp_polynomial_t* A) {
  return coefficient_is_univariate(&A->data);
}

int lp_polynomial_is_univariate_m(const lp_polynomial_t* A, const lp_assignment_t* m) {
  if (lp_polynomial_is_constant(A)) {
    return 0;
  }
  lp_variable_t top = lp_polynomial_top_variable(A);
  if (lp_assignment_get_value(m, top)->type != LP_VALUE_NONE) {
    return 0;
  }
  lp_variable_list_t vars;
  lp_variable_list_construct(&vars);
  lp_polynomial_get_variables(A, &vars);
  size_t i;
  for (i = 0; i < vars.list_size; ++ i) {
    lp_variable_t x = vars.list[i];
    if (x != top && lp_assignment_get_value(m, x)->type == LP_VALUE_NONE) {
      break;
    }
  }
  if (i < vars.list_size) {
    return 0;
  }
  lp_variable_list_destruct(&vars);
  return 1;
}

lp_upolynomial_t* lp_polynomial_to_univariate(const lp_polynomial_t* A) {
  if (!coefficient_is_univariate(&A->data)) {
    return 0;
  } else {
    return coefficient_to_univariate(A->ctx, &A->data);
  }
}

int lp_polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);

  // Check that all variables are assigned
  lp_variable_list_t A_vars;
  lp_variable_list_construct(&A_vars);
  coefficient_get_variables(&A->data, &A_vars);
  size_t i;
  for (i = 0; i < A_vars.list_size; ++ i) {
    lp_variable_t x = A_vars.list[i];
    if (lp_assignment_get_value(m, x)->type == LP_VALUE_NONE) {
      return -2;
    }
  }
  lp_variable_list_destruct(&A_vars);

  return coefficient_sgn(A->ctx, &A->data, m);
}

lp_value_t* lp_polynomial_evaluate(const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);

  // Check that all variables are assigned
  lp_variable_list_t A_vars;
  lp_variable_list_construct(&A_vars);
  coefficient_get_variables(&A->data, &A_vars);
  size_t i;
  for (i = 0; i < A_vars.list_size; ++ i) {
    lp_variable_t x = A_vars.list[i];
    if (lp_assignment_get_value(m, x)->type == LP_VALUE_NONE) {
      return 0;
    }
  }
  lp_variable_list_destruct(&A_vars);

  return coefficient_evaluate(A->ctx, &A->data, m);
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

void lp_polynomial_add_monomial(lp_polynomial_t* S, const lp_monomial_t* M) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_add("); lp_polynomial_print(S, trace_out); tracef(", "); lp_monomial_print(S->ctx, M, trace_out); tracef(")\n");
    lp_variable_order_print(
        S->ctx->var_order, S->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(S);

  coefficient_add_monomial(S->ctx, &S->data, M);

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

/** Compute P = A1 * C. */
void lp_polynomial_mul_integer(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_integer_t* C) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mul_c("); lp_polynomial_print(P, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_integer_print(C, trace_out); tracef(")\n");
    lp_variable_order_print(
        A1->ctx->var_order, A1->ctx->var_db,
        trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A1);

  lp_polynomial_set_context(P, A1->ctx);

  coefficient_mul_integer(P->ctx, &P->data, &A1->data, C);

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

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("A = "); lp_polynomial_print(A, trace_out); tracef("\n");
    tracef("B = "); lp_polynomial_print(B, trace_out); tracef("\n");
    tracef("var = %s\n", lp_variable_db_get_name(A->ctx->var_db, lp_polynomial_top_variable(A)));
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out); tracef("\n");
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

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("Sq Factor A = "); lp_polynomial_print(A, trace_out); tracef("\n");
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out); tracef("\n");
  }

//  assert(*factors == 0);
//  assert(*multiplicities == 0);
//  assert(*size == 0);

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

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("Sq Factor: result size = %zu\n", *size);
  }

  coefficient_factors_destruct(&coeff_factors);
}

void lp_polynomial_roots_isolate(const lp_polynomial_t* A, const lp_assignment_t* M, lp_value_t* roots, size_t* roots_size) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_roots_isolate("); lp_polynomial_print(A, trace_out); tracef(")\n");
  }

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("A = "); lp_polynomial_print(A, trace_out); tracef("\n");
    tracef("var = %s\n", lp_variable_db_get_name(A->ctx->var_db, lp_polynomial_top_variable(A)));
    lp_assignment_print(M, trace_out); tracef("\n");
  }


  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  lp_variable_t x = lp_polynomial_top_variable(A);
  assert(x != lp_variable_null);

  lp_value_t x_value_backup;
  if (lp_assignment_get_value(M, x)->type != LP_VALUE_NONE) {
    lp_value_construct_copy(&x_value_backup, lp_assignment_get_value(M, x));
    lp_assignment_set_value((lp_assignment_t*) M, x, 0);
  } else {
    lp_value_construct_none(&x_value_backup);
  }

  size_t i;

  lp_polynomial_t** factors = 0;
  size_t* multiplicities = 0;
  size_t factors_size = 0;

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_roots_isolate(): factoring\n");
  }

  // Get the square-free factorization
  lp_polynomial_factor_square_free(A, &factors, &multiplicities, &factors_size);

  // Count the max number of roots
  size_t total_degree = 0;
  size_t factor_i;
  for (factor_i = 0; factor_i < factors_size; ++ factor_i) {
    // Add the degree of the polynomial if in top variable
    if (lp_polynomial_top_variable(factors[factor_i]) == x) {
      total_degree += lp_polynomial_degree(factors[factor_i]);
    }
  }

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_roots_isolate(): factors = %zu, total_degree = %zi\n", factors_size, total_degree);
  }

  // Allocate enough space for the roots
  lp_value_t* roots_tmp = malloc(sizeof(lp_value_t)*total_degree);
  size_t roots_tmp_size = 0;

  for (factor_i = 0; factor_i < factors_size; ++ factor_i) {
    // The factor we are working with
    const lp_polynomial_t* factor = factors[factor_i];
    // Get the roots if not a constant
    if (x == lp_polynomial_top_variable(factor)) {
      // Proper polynomial in x
      assert(roots_tmp_size + lp_polynomial_degree(factor) <= total_degree);
      lp_value_t* current_roots = roots_tmp + roots_tmp_size;
      size_t current_roots_size = 0;
      coefficient_roots_isolate(A->ctx, &factor->data, M, current_roots, &current_roots_size);
      roots_tmp_size += current_roots_size;
      assert(roots_tmp_size <= total_degree);
    } else {
      // Polynomial in some other variable -- we need to check the sign: if 0
      // then there is no roots all together
      int sgn = lp_polynomial_sgn(factor, M);
      if (sgn == 0) {
        for (i = 0; i < roots_tmp_size; ++ i) {
          lp_value_destruct(roots_tmp + i);
        }
        roots_tmp_size = 0;
        break;
      }
    }
  }

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_root_isolate("); lp_polynomial_print(A, trace_out); tracef("): unsorted roots\n")
    for (i = 0; i < roots_tmp_size; ++ i) {
      tracef("%zu :", i); lp_value_print(roots_tmp + i, trace_out); tracef("\n");
    }
  }

  if (roots_tmp_size > 0) {
    // Sort the roots
    qsort(roots_tmp, roots_tmp_size, sizeof(lp_value_t), lp_value_cmp_void);

    if (trace_is_enabled("polynomial")) {
      tracef("polynomial_root_isolate("); lp_polynomial_print(A, trace_out); tracef("): sorted roots\n")
      for (i = 0; i < roots_tmp_size; ++ i) {
        tracef("%zu :", i); lp_value_print(roots_tmp + i, trace_out); tracef("\n");
      }
    }

    // Remove any duplicates
    size_t to_keep;
    for (to_keep = 1, i = 1; i < roots_tmp_size; ++ i) {
      if (lp_value_cmp(roots_tmp + i, roots_tmp  + to_keep-1) != 0) {
        // If different copy over
        if (i != to_keep) {
          lp_value_assign(roots_tmp + to_keep, roots_tmp + i);
        }
        // This one is a keeper
        to_keep ++;
      }
    }
    for (i = to_keep; i < roots_tmp_size; ++ i) {
      lp_value_destruct(roots_tmp + i);
    }
    roots_tmp_size = to_keep;

    // Copy over the roots
    memcpy(roots, roots_tmp, roots_tmp_size*sizeof(lp_value_t));
  }

  // Set the new size
  *roots_size = roots_tmp_size;

  // Reset the value
  if (x_value_backup.type != LP_VALUE_NONE) {
    lp_assignment_set_value((lp_assignment_t*) M, x, &x_value_backup);
  }

  // Destroy the temps
  for (i = 0; i < factors_size; ++ i) {
    lp_polynomial_destruct(factors[i]);
    free(factors[i]);
  }
  free(factors);
  free(multiplicities);
  free(roots_tmp);
  lp_value_destruct(&x_value_backup);
}

lp_feasibility_set_t* lp_polynomial_constraint_get_feasible_set(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, int negated, const lp_assignment_t* M) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_get_feasible_set("); lp_polynomial_print(A, trace_out); tracef(", "); lp_sign_condition_print(sgn_condition, trace_out); tracef(")\n");
  }

  assert(!lp_polynomial_is_constant(A));

  // Make sure we're in the right order
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  // Make sure that the top variable is unassigned
  lp_variable_t x = coefficient_top_variable(&A->data);
  assert(lp_assignment_get_value(M, x)->type == LP_VALUE_NONE);

  // Get the degree of the polynomial, respecting the model
  size_t degree = coefficient_degree_m(A->ctx, &A->data, M);

  // Negate the constraint if negated
  if (negated) {
    sgn_condition = lp_sign_condition_negate(sgn_condition);
  }

  if (degree == 0) {
    // Evaluates to constant
    int sgn = coefficient_sgn(A->ctx, coefficient_get_coefficient_safe(A->ctx, &A->data, 0, x), M);

    if (trace_is_enabled("polynomial")) {
      tracef("polynomial_get_feasible_set(");
      lp_polynomial_print(A, trace_out);
      tracef(", "); lp_sign_condition_print(sgn_condition, trace_out);
      tracef(") => evaluates to constant of sign %d\n", sgn);
    }

    if (lp_sign_condition_consistent(sgn_condition, sgn)) {
      // Consistent for any x
      return lp_feasibility_set_new_full();
    } else {
      // No x
      return lp_feasibility_set_new_internal(0);
    }
  }

  // Get the roots of the polynomial
  size_t roots_size;
  lp_value_t* roots = malloc(sizeof(lp_value_t)*degree);
  lp_polynomial_roots_isolate(A, M, roots, &roots_size);

  //
  // We have roots:
  //
  //   r0   r1   ...   r_n
  //
  // We start at -inf and go to -inf, where we can compute the sign of the
  // polynomial.
  //

  size_t signs_size = 2*roots_size + 1;
  int* signs = malloc(sizeof(int)*signs_size);

  // Get the first non-vanishing coefficient, or constant otherwise
  int sgn_lc = coefficient_sgn(A->ctx, coefficient_get_coefficient(&A->data, degree), M);

  // Signs at -inf and +inf
  signs[0] = degree % 2 ? -sgn_lc : sgn_lc;
  signs[signs_size-1] = sgn_lc;

  // Signs at root
  size_t i;
  lp_value_t m;
  lp_value_construct_none(&m);
  for (i = 0; i < roots_size; ++ i) {
    signs[2*i+1] = 0;
    if (i+1<roots_size) {
      lp_value_get_value_between(roots + i, 1, roots + i + 1, 1, &m);
      lp_assignment_set_value((lp_assignment_t*) M, x, &m);
      signs[2*i+2] = coefficient_sgn(A->ctx, &A->data, M);
      lp_assignment_set_value((lp_assignment_t*) M, x, 0);
    }
  }
  lp_value_destruct(&m);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_get_feasible_set():");
    for (i = 0; i < signs_size; ++ i) {
      if (signs[i] < 0) { tracef(" -"); }
      else if (signs[i] == 0) { tracef(" 0"); }
      else { tracef(" +"); }
    }
    tracef("\n");
  }

  // Count the number of intervals
  size_t intervals_size = 0, lb, ub;
  for (lb = 0; lb < signs_size; ) {
    // Find lower bound
    for (; lb < signs_size && !lp_sign_condition_consistent(sgn_condition, signs[lb]); lb ++) {}
    if (lb < signs_size) {
      // Found one
      intervals_size ++;
      // Find the upper bound
      for (ub = lb + 1; ub < signs_size && lp_sign_condition_consistent(sgn_condition, signs[ub]); ub ++) {}
      // Continue with the next one
      lb = ub;
    }
  }


  lp_feasibility_set_t* result = lp_feasibility_set_new_internal(intervals_size);

  lp_value_t inf_neg, inf_pos;
  lp_value_construct(&inf_neg, LP_VALUE_MINUS_INFINITY, 0);
  lp_value_construct(&inf_pos, LP_VALUE_PLUS_INFINITY, 0);

  // Go through signs and collect the contiguous intervals
  size_t interval = 0;
  for (lb = 0; lb < signs_size; ) {
    // find lower bound
    for (; lb < signs_size && !lp_sign_condition_consistent(sgn_condition, signs[lb]); lb ++) {}
    if (lb < signs_size) {
      // find upper bound
      for (ub = lb + 1; ub < signs_size && lp_sign_condition_consistent(sgn_condition, signs[ub]); ub ++) {}

      // Found the interval
      if (lb == (ub + 1) && lb % 2) {
        // it's a point
        lp_interval_construct_point(result->intervals + interval, roots + (lb / 2));
      } else {
        // It's an interval
        const lp_value_t* lb_value;
        const lp_value_t* ub_value;
        int lb_strict, ub_strict;

        // lb point to the first satisfied
        // ub is last satisfied +1

        // Lower bound
        if (lb % 2 == 1) {
          // Root
          lb_value = roots + (lb/2);
          lb_strict = 0;
        } else {
          // Interval
          if (lb == 0) {
            // -inf
            lb_value = &inf_neg;
          } else {
            // Root bounding the interval
            lb_value = roots + ((lb-1)/2);
          }
          lb_strict = 1;
        }

        // Upper bound
        if (ub % 2 == 0) {
          // Root
          ub_value = roots + ((ub-1)/2);
          ub_strict = 0;
        } else {
          // Interval
          if (ub == signs_size) {
            // +inf
            ub_value = &inf_pos;
          } else {
            // Root bounding the interval
            ub_value = roots + (ub/2);
          }
          ub_strict = 1;
        }

        // Construct the interval
        lp_interval_construct(result->intervals + interval, lb_value, lb_strict, ub_value, ub_strict);
      }

      // Done with this interval, continue with the next
      interval ++;
      lb = ub;
    }
  }

  assert(interval == intervals_size);
  assert(interval == result->capacity);
  result->size = intervals_size;

  // Remove the signs array
  free(signs);

  // Remove the roots
  for (i = 0; i < roots_size; ++ i) {
    lp_value_destruct(roots + i);
  }
  free(roots);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_get_feasible_set(");
    lp_polynomial_print(A, trace_out);
    tracef(", "); lp_sign_condition_print(sgn_condition, trace_out);
    tracef(") => "); lp_feasibility_set_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_feasibility_set_t* lp_polynomial_root_constraint_get_feasible_set(const lp_polynomial_t* A, size_t root_index, lp_sign_condition_t sgn_condition, int negated, const lp_assignment_t* M) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_get_feasible_set_root("); lp_polynomial_print(A, trace_out); tracef(", %zu, ", root_index); lp_sign_condition_print(sgn_condition, trace_out); tracef(")\n");
  }

  assert(!lp_polynomial_is_constant(A));

  // Make sure we're in the right order
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  // Make sure that the top variable is unassigned
  assert(coefficient_top_variable(&A->data) != lp_variable_null);

  // Get the degree of the polynomial, respecting the model
  size_t degree = coefficient_degree_m(A->ctx, &A->data, M);
  if (degree == 0) {
    // Polynomial evaluates to 0, no roots => root_index < 0 is false
    if (!negated) {
      return lp_feasibility_set_new_internal(0);
    } else {
      return lp_feasibility_set_new_full();
    }
  }

  // Get the roots of the polynomial
  size_t roots_size;
  lp_value_t* roots = malloc(sizeof(lp_value_t)*degree);
  lp_polynomial_roots_isolate(A, M, roots, &roots_size);

  /**
   * Example:
   *
   *   x <_r root(k, p(x, y))
   *
   * not negated:
   *
   *   k < rootcount(x, p(x, y)) && x < root(k, p(x, y))
   *
   *   if (k < rootcount): x \in (-inf, root(k, p(x, y)))
   *   else              : x \in {}
   *
   * negated:
   *
   *   k >= rootcount(x, p(x, y)) || x >= root(k, p(x, y))
   *
   *   if (k >= rootcount): x in (-inf, +inf)
   *   else               : x in (root(k, p(x, y), +inf)
   *
   */

  lp_feasibility_set_t* result = 0;
  if (root_index >= roots_size) {
    // Just empty
    if (!negated) {
      result = lp_feasibility_set_new_internal(0);
    } else {
      result = lp_feasibility_set_new_full();
    }
  } else {

    if (negated) {
      sgn_condition = lp_sign_condition_negate(sgn_condition);
    }

    lp_value_t inf_pos, inf_neg;
    lp_value_construct(&inf_pos, LP_VALUE_PLUS_INFINITY, 0);
    lp_value_construct(&inf_neg, LP_VALUE_MINUS_INFINITY, 0);

    switch (sgn_condition) {
    case LP_SGN_LT_0:
      // (-inf, root)
      result = lp_feasibility_set_new_internal(1);
      lp_interval_construct(result->intervals, &inf_neg, 1, roots + root_index, 1);
      result->size = 1;
      break;
    case LP_SGN_LE_0:
      // (-inf, root]
      result = lp_feasibility_set_new_internal(1);
      lp_interval_construct(result->intervals, &inf_neg, 1, roots + root_index, 0);
      result->size = 1;
      break;
    case LP_SGN_EQ_0:
      // [root, root]
      result = lp_feasibility_set_new_internal(1);
      lp_interval_construct_point(result->intervals, roots + root_index);
      result->size = 1;
      break;
    case LP_SGN_NE_0:
      // (-inf, root) (root, +inf)
      result = lp_feasibility_set_new_internal(2);
      lp_interval_construct(result->intervals, &inf_neg, 1, roots + root_index, 1);
      lp_interval_construct(result->intervals + 1, roots + root_index, 1, &inf_pos, 1);
      result->size = 2;
      break;
    case LP_SGN_GT_0:
      // (root, +inf)
      result = lp_feasibility_set_new_internal(1);
      lp_interval_construct(result->intervals, roots + root_index, 1, &inf_pos, 1);
      result->size = 1;
      break;
    case LP_SGN_GE_0:
      // [root, +inf)
      result = lp_feasibility_set_new_internal(1);
      lp_interval_construct(result->intervals, roots + root_index, 0, &inf_pos, 1);
      result->size = 1;
      break;
    }

    lp_value_destruct(&inf_neg);
    lp_value_destruct(&inf_pos);
  }

  // Remove the roots
  size_t i;
  for (i = 0; i < roots_size; ++ i) {
    lp_value_destruct(roots + i);
  }
  free(roots);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_get_feasible_set(");
    lp_polynomial_print(A, trace_out);
    tracef(", "); lp_sign_condition_print(sgn_condition, trace_out);
    tracef(") => "); lp_feasibility_set_print(result, trace_out); tracef("\n");
  }

  return result;

}

void lp_polynomial_get_variables(const lp_polynomial_t* A, lp_variable_list_t* vars) {
  coefficient_get_variables(&A->data, vars);
}

size_t lp_polynomial_hash(const lp_polynomial_t* A_const) {
  if (!A_const->hash) {
    size_t hash = coefficient_hash(A_const->ctx, &A_const->data);
    if (hash == 0) {
      hash ++;
    }
    ((lp_polynomial_t*)A_const)->hash = hash;
  }
  return A_const->hash;
}

int lp_polynomial_eq(const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

  size_t A1_hash = lp_polynomial_hash(A1);
  size_t A2_hash = lp_polynomial_hash(A2);

  if (A1_hash != A2_hash) {
    // Different hashes => different
    return 0;
  }

  // Compare
  return lp_polynomial_cmp(A1, A2) == 0;
}

void lp_polynomial_traverse(const lp_polynomial_t* A, lp_polynomial_traverse_f f, void* data) {

  lp_polynomial_external_clean(A);

  lp_monomial_t m;
  lp_monomial_construct(A->ctx, &m);
  coefficient_traverse(A->ctx, &A->data, f, &m, data);
  lp_monomial_destruct(&m);
}

int lp_polynomial_constraint_evaluate(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, const lp_assignment_t* M) {

  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  // Evaluate the sign and check
  int p_sign = lp_polynomial_sgn(A, M);
  return lp_sign_condition_consistent(sgn_condition, p_sign);
}

int lp_polynomial_root_constraint_evaluate(const lp_polynomial_t* A, size_t root_index, lp_sign_condition_t sgn_condition, const lp_assignment_t* M) {


  int eval;

  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  lp_variable_t x = lp_polynomial_top_variable(A);
  assert(x != lp_variable_null);

  size_t degree = lp_polynomial_degree(A);
  lp_value_t* roots = malloc(sizeof(lp_value_t)*degree);
  size_t roots_size = 0;

  // Get the root we're interested in
  lp_polynomial_roots_isolate(A, M, roots, &roots_size);
  if (root_index < roots_size) {
    // Compare
    const lp_value_t* x_value = lp_assignment_get_value(M, x);
    int cmp = lp_value_cmp(x_value, roots + root_index);
    eval = lp_sign_condition_consistent(sgn_condition, cmp);
  } else {
    // false, not enough roots
    eval = 0;
  }

  size_t i;
  for (i = 0; i < roots_size; ++ i) {
    lp_value_destruct(roots + i);
  }
  free(roots);

  return eval;
}
