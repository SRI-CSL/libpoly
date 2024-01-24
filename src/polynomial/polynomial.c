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

#include <upolynomial.h>
#include <feasibility_set.h>
#include <variable_db.h>
#include <variable_list.h>

#include "polynomial/polynomial.h"
#include "polynomial/coefficient.h"

#include "polynomial/gcd.h"
#include "polynomial/factorization.h"
#include "polynomial/output.h"

#include "number/rational.h"
#include "number/integer.h"

#include "polynomial/feasibility_set.h"
#include "polynomial/polynomial_vector.h"

#include "utils/debug_trace.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define SWAP(type, x, y) ({ type tmp = x; x = y; y = tmp; })

static
void check_polynomial_assignment(const lp_polynomial_t* A, const lp_assignment_t* M, lp_variable_t x) {

  // Given variable must be top
  if (A->data.type != COEFFICIENT_NUMERIC && x != lp_variable_null && VAR(&A->data) != x) {
    fprintf(stderr, "M = "); lp_assignment_print(M, stderr); fprintf(stderr, "\n");
    fprintf(stderr, "A = "); lp_polynomial_print(A, stderr); fprintf(stderr, "\n");
    assert(0);
  }

  lp_variable_list_t vars;
  lp_variable_list_construct(&vars);
  lp_polynomial_get_variables(A, &vars);
  size_t i = 0;
  for (i = 0; i < vars.list_size; ++i) {
    lp_variable_t y = vars.list[i];
    if (x != y && lp_assignment_get_value(M, y)->type == LP_VALUE_NONE) {
      lp_polynomial_print(A, trace_out);
      tracef("\n");
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

int lp_polynomial_check_order(const lp_polynomial_t* A) {
  return coefficient_in_order(A->ctx, &A->data);
}

void lp_polynomial_ensure_order(lp_polynomial_t* A) {
  coefficient_order(A->ctx, &A->data);
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

lp_polynomial_t* lp_polynomial_new_from_coefficient(const lp_polynomial_context_t* ctx, const coefficient_t* from) {
  lp_polynomial_t* result = lp_polynomial_alloc();
  lp_polynomial_construct_from_coefficient(result, ctx, from);
  return result;
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

int lp_polynomial_is_linear(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_is_linear(&A->data);
}

int lp_polynomial_lc_is_constant(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_lc(&A->data)->type == COEFFICIENT_NUMERIC;
}

int lp_polynomial_lc_sgn(const lp_polynomial_t* A) {
  lp_polynomial_external_clean(A);
  return coefficient_lc_sgn(A->ctx, &A->data);
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
  coefficient_reductum_m(A->ctx, &R->data, &A->data, m, 0);
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

  int result = (i == vars.list_size);
  lp_variable_list_destruct(&vars);
  return result;
}

lp_upolynomial_t* lp_polynomial_to_univariate(const lp_polynomial_t* A) {
  if (!(coefficient_is_constant(&A->data) || coefficient_is_univariate(&A->data))) {
    return NULL;
  } else {
    return coefficient_to_univariate(A->ctx, &A->data);
  }
}

lp_upolynomial_t* lp_polynomial_to_univariate_m(const lp_polynomial_t* A, const lp_assignment_t* m) {
  return coefficient_to_univariate_m(A->ctx, &A->data, m);
}

int lp_polynomial_is_assigned(const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);
  return coefficient_is_assigned(A->ctx, &A->data, m);
}

int lp_polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, m, lp_variable_null);
  }

  return coefficient_sgn(A->ctx, &A->data, m);
}

void lp_polynomial_interval_value(const lp_polynomial_t* A, const lp_interval_assignment_t* m, lp_interval_t* result) {
  lp_polynomial_external_clean(A);
  coefficient_interval_value(A->ctx, &A->data, m, result);
}

lp_value_t* lp_polynomial_evaluate(const lp_polynomial_t* A, const lp_assignment_t* m) {
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, m, lp_variable_null);
  }

  return coefficient_evaluate(A->ctx, &A->data, m);
}

void lp_polynomial_evaluate_integer(const lp_polynomial_t* A, const lp_assignment_t* m, lp_integer_t *value) {
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, m, lp_variable_null);
  }

  coefficient_evaluate_integer(A->ctx, &A->data, m, value);
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

void lp_polynomial_pdivrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

    if (trace_is_enabled("polynomial")) {
        tracef("polynomial_pdirvrem("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
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

    coefficient_pdivrem(D->ctx, &D->data, &R->data, &A1->data, &A2->data);

    if (trace_is_enabled("polynomial")) {
        tracef("polynomial_pdirvrem() => ("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(")\n");
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

void lp_polynomial_spdivrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2) {

    if (trace_is_enabled("polynomial")) {
        tracef("lp_polynomial_spdirvrem("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(", "); lp_polynomial_print(A1, trace_out); tracef(", "); lp_polynomial_print(A2, trace_out); tracef(")\n");
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

    coefficient_spdivrem(D->ctx, &D->data, &R->data, &A1->data, &A2->data);

    if (trace_is_enabled("polynomial")) {
        tracef("lp_polynomial_spdirvrem() => ("); lp_polynomial_print(D, trace_out); tracef(", "); lp_polynomial_print(R, trace_out); tracef(")\n");
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

  coefficient_divrem(D->ctx, &D->data, &R->data, &A1->data, &A2->data);

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

  coefficient_reduce(ctx, &A->data, &B->data, &P->data, &Q->data, &R->data, REMAINDERING_PSEUDO_DENSE);

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_reduce() =>\n");
    tracef("\t P = "); lp_polynomial_print(P, trace_out); tracef("\n");
    tracef("\t Q = "); lp_polynomial_print(Q, trace_out); tracef("\n");
    tracef("\t R = "); lp_polynomial_print(R, trace_out); tracef("\n");
  }
}

void lp_polynomial_cont(lp_polynomial_t* cont, const lp_polynomial_t* A) {
  const lp_polynomial_context_t* ctx = A->ctx;
  lp_polynomial_external_clean(A);
  lp_polynomial_set_context(cont, ctx);
  coefficient_cont(ctx, &cont->data, &A->data);
}

void lp_polynomial_pp(lp_polynomial_t* pp, const lp_polynomial_t* A) {
  const lp_polynomial_context_t* ctx = A->ctx;
  lp_polynomial_external_clean(A);
  lp_polynomial_set_context(pp, ctx);
  coefficient_pp(ctx, &pp->data, &A->data);
}

void lp_polynomial_pp_cont(lp_polynomial_t* pp, lp_polynomial_t* cont, const lp_polynomial_t* A) {
  const lp_polynomial_context_t* ctx = A->ctx;
  lp_polynomial_external_clean(A);
  lp_polynomial_set_context(pp, ctx);
  lp_polynomial_set_context(cont, ctx);
  coefficient_pp_cont(ctx, &pp->data, &cont->data, &A->data);
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
  // the last one is not ignored right now
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

void lp_polynomial_srs(lp_polynomial_t** srs, const lp_polynomial_t* A, const lp_polynomial_t* B) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_srs("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(")\n");
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
    lp_polynomial_srs(srs, B, A);
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
  coefficient_t*srs_coeff = malloc(sizeof(coefficient_t)*size);
  size_t i;
  for (i = 0; i < size; ++ i) {
    coefficient_construct(ctx, srs_coeff + i);
  }

  // Compute
  coefficient_srs(ctx, srs_coeff, &A->data, &B->data);

  // Construct the output (one less, we ignore the final 1)
  // the last one is not ignored right now
  for (i = 0; i < size; ++ i) {
    lp_polynomial_t tmp;
    lp_polynomial_construct_from_coefficient(&tmp, ctx, srs_coeff + i);
    lp_polynomial_swap(&tmp, srs[i]);
    lp_polynomial_destruct(&tmp);
    coefficient_destruct(&srs_coeff[i]);
  }

  free(srs_coeff);

  if (trace_is_enabled("polynomial")) {
    for (i = 0; i < size; ++ i) {
      tracef("SRS[%zu] = ", i); lp_polynomial_print(srs[i], trace_out); tracef("\n");
    }
  }
}

lp_polynomial_vector_t* lp_polynomial_mgcd(const lp_polynomial_t* A, const lp_polynomial_t* B, const lp_assignment_t* m) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_mgcd("); lp_polynomial_print(A, trace_out); tracef(", "); lp_polynomial_print(B, trace_out); tracef(")\n");
  }

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  const lp_polynomial_context_t* ctx = A->ctx;
  assert(lp_polynomial_context_equal(B->ctx, ctx));

  lp_polynomial_external_clean(A);
  lp_polynomial_external_clean(B);

  // Compute it
  return coefficient_mgcd(ctx, &A->data, &B->data, m);
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

void lp_polynomial_factor_content_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_factor_content_free("); lp_polynomial_print(A, trace_out); tracef(")\n");
  }

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("Content Factor A = "); lp_polynomial_print(A, trace_out); tracef("\n");
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out); tracef("\n");
  }

  const lp_polynomial_context_t* ctx = A->ctx;

  if (trace_is_enabled("polynomial")) {
    lp_variable_order_print(A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  lp_polynomial_external_clean(A);
  coefficient_factors_t coeff_factors;
  coefficient_factors_construct(&coeff_factors);

  coefficient_factor_content_free(ctx, &A->data, &coeff_factors);

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
    tracef("Content Factor: result size = %zu\n", *size);
  }

  coefficient_factors_destruct(&coeff_factors);
}


void lp_polynomial_roots_isolate(const lp_polynomial_t* A, const lp_assignment_t* M, lp_value_t* roots, size_t* roots_size) {

  if (trace_is_enabled("polynomial")) {
    tracef("polynomial_roots_isolate("); lp_polynomial_print(A, trace_out); tracef(")\n");
  }

  if (trace_is_enabled("polynomial::expensive")) {
    tracef("Root isolation\n");
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

  // Get the reduced polynomial
  lp_polynomial_t A_r;
  lp_polynomial_construct(&A_r, A->ctx);
  lp_polynomial_reductum_m(&A_r, A, M);
  assert(x == lp_polynomial_top_variable(A));

  // Get the square-free factorization
  lp_polynomial_factor_square_free(&A_r, &factors, &multiplicities, &factors_size);

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
    tracef("polynomial_root_isolate("); lp_polynomial_print(A, trace_out); tracef("): unsorted roots\n");
    for (i = 0; i < roots_tmp_size; ++ i) {
      tracef("%zu :", i); lp_value_print(roots_tmp + i, trace_out); tracef("\n");
    }
  }

  if (roots_tmp_size > 0) {
    // Sort the roots
    qsort(roots_tmp, roots_tmp_size, sizeof(lp_value_t), lp_value_cmp_void);

    if (trace_is_enabled("polynomial")) {
      tracef("polynomial_root_isolate("); lp_polynomial_print(A, trace_out); tracef("): sorted roots\n");
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
  lp_polynomial_destruct(&A_r);
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

  // Top variable
  lp_variable_t x = coefficient_top_variable(&A->data);

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

int lp_polynomial_constraint_infer_bounds(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, int negated, lp_interval_assignment_t* M) {

  // Negate the constraint if negated
  if (negated) {
    sgn_condition = lp_sign_condition_negate(sgn_condition);
  }

  const lp_polynomial_context_t* ctx = A->ctx;

  switch (sgn_condition) {
  case LP_SGN_LT_0: // |x| - d < 0 => d < x < d
  case LP_SGN_LE_0: // |x| - d <= 0 => d <= x <= d
    break;
  case LP_SGN_EQ_0: {
    // |x| - d == 0 == both <=, >=
    int r = lp_polynomial_constraint_infer_bounds(A, LP_SGN_LE_0, 0, M);
    if (r) {
      return r;
    }
    return lp_polynomial_constraint_infer_bounds(A, LP_SGN_LE_0, 0, M);
  }
  case LP_SGN_NE_0: // |x| - d != 0 => ?
    return 0;
  case LP_SGN_GT_0: {
    lp_polynomial_t* A_neg = lp_polynomial_new(ctx);
    lp_polynomial_neg(A_neg, A);
    int result = lp_polynomial_constraint_infer_bounds(A_neg, LP_SGN_LT_0, 0, M);
    lp_polynomial_delete(A_neg);
    return result;
  }
  case LP_SGN_GE_0: {
    lp_polynomial_t* A_neg = lp_polynomial_new(ctx);
    lp_polynomial_neg(A_neg, A);
    int result = lp_polynomial_constraint_infer_bounds(A_neg, LP_SGN_LE_0, 0, M);
    lp_polynomial_delete(A_neg);
    return result;
  }
  }

  if (trace_is_enabled("polynomial::bounds")) {
    tracef("lp_polynomial_constraint_infer_bounds("); lp_polynomial_print(A, trace_out); tracef(", "); lp_sign_condition_print(sgn_condition, trace_out); tracef(")\n");
  }

  // Make sure we're in the right order
  lp_polynomial_external_clean(A);

  // For each variable in A, see A(x). If
  //
  //    A(x) = Ax^2 + Bx + C,
  //
  // with A and B appropriate constants, then we can write this as
  //
  //   A(x) = a^2 x - 2ab x + b^2 + C - b^2
  //        = (ax - b)^2 + C - b^2
  //
  // If this is the case, and we can recursively rewrite the C further into
  //
  //  (a1 x1 - b1)^2 + ... + (an xn - bn)^2 + C0 - b1^2 - ... - bn^2
  //
  // With D = (b1^2 + ... + bn^2) - C0, and d = sqrt(D), we can infer bounds
  //
  //    (ak x - bk)^2 <= D => -(d + bk)/ak <= x <= (d + bk)/ak
  //
  // In fact, we can just solve
  //
  //    ak^2 - 2akabk x - D^2 == 0
  //
  // and make the interval (r1, r2), [r1], or 0, depending on the number of
  // roots and the sign condition.
  //
  // We can do this easily in the original polynomial, but D is rational so we need
  // to multiply it out.

  lp_integer_t tmp_z, B_sq, A4;
  integer_construct(&tmp_z);
  integer_construct(&B_sq);
  integer_construct(&A4);

  lp_rational_t D;
  rational_construct(&D);

  int ok = 1;
  const coefficient_t* Ak = &A->data;
  while (ok && Ak->type != COEFFICIENT_NUMERIC) {
    // Check if Ak = (ax - b)^2 + Ak-1 = Ax^2 - Bx + ...
    if (trace_is_enabled("polynomial::bounds")) {
      tracef("A_k = "); coefficient_print(ctx, Ak, trace_out); tracef("\n");
      tracef("D = "); rational_print(&D, trace_out); tracef("\n");
    }
    size_t Ak_degree = coefficient_degree(Ak);
    if (Ak_degree == 2) {
      const coefficient_t* A = coefficient_get_coefficient(Ak, 2);
      const coefficient_t* B = coefficient_get_coefficient(Ak, 1);
      if (A->type != COEFFICIENT_NUMERIC || B->type != COEFFICIENT_NUMERIC) {
        ok = 0;
      } else if (integer_sgn(lp_Z, &A->value.num) > 0) {
        //  A(x) = (ax - b)^2 + ...
        //  A = a^2, B=-2ab
        //  A(x) = Ax^2 + Bx + (b^2 + ...)
        //  b = -B/2*a
        //  b^2 = B^2/4*A
        integer_mul(lp_Z, &B_sq, &B->value.num, &B->value.num);
        integer_mul_int(lp_Z, &A4, &A->value.num, 4);
        lp_rational_t tmp_q;
        rational_construct_from_div(&tmp_q, &B_sq, &A4);
        rational_add(&D, &D, &tmp_q);
        rational_destruct(&tmp_q);
        // Go to the next one
        Ak = COEFF(Ak, 0);
      } else {
        // Cannot do square root, must be positive
        ok = 0;
      }
    } else {
      ok = 0;
    }
  }

  int conflict = 0;
  if (ok) {
    // D = d^2
    lp_rational_t tmp_q;
    rational_construct_from_integer(&tmp_q, &Ak->value.num);
    rational_sub(&D, &D, &tmp_q);
    rational_destruct(&tmp_q);

    // Construct a polynomial for root finding
    coefficient_t f, next;
    coefficient_construct_copy(ctx, &next, &A->data);
    coefficient_construct(ctx, &f);

    // Now, we traverse again, and get the intervals
    while (next.type != COEFFICIENT_NUMERIC) {
      // Continue to the next one
      coefficient_swap(&f, &next);
      coefficient_assign_int(ctx, &next, 0);
      coefficient_swap(&next, COEFF(&f, 0));

      if (trace_is_enabled("polynomial::bounds")) {
        tracef("f = "); coefficient_print(ctx, &f, trace_out); tracef("\n");
        tracef("D = "); rational_print(&D, trace_out); tracef("\n");
      }
      // Solve Ax^2 + Bx + b^2 - D == 0
      // D is rational p/q we solve
      // B = -2ab => b^2 = B^2/4a^2 = B^2/4A
      // b^2 = B^2/4*A
      const coefficient_t* A = COEFF(&f, 2);
      const coefficient_t* B = COEFF(&f, 1);
      integer_mul(lp_Z, &B_sq, &B->value.num, &B->value.num);
      integer_mul_int(lp_Z, &A4, &A->value.num, 4);
      rational_construct_from_div(&tmp_q, &B_sq, &A4);
      rational_sub(&tmp_q, &tmp_q, &D);
      // Add p/q to polynomial to solve
      const lp_integer_t* p = rational_get_num_ref(&tmp_q);
      const lp_integer_t* q = rational_get_den_ref(&tmp_q);
      coefficient_assign_int(ctx, COEFF(&f, 0), 0);
      coefficient_mul_integer(ctx, &f, &f, q);
      coefficient_assign_integer(ctx, COEFF(&f, 0), p);
      rational_destruct(&tmp_q);
      if (trace_is_enabled("polynomial::bounds")) {
        tracef("f = "); coefficient_print(ctx, &f, trace_out); tracef("\n");
      }
      // Get the roots
      lp_value_t roots[2];
      size_t roots_size = 0;
      coefficient_roots_isolate_univariate(ctx, &f, roots, &roots_size);
      // Create the result
      lp_variable_t x = VAR(&f);
      lp_interval_t x_interval;
      if (roots_size == 0) {
        // No roots, inconsistent
        conflict = 1;
        break;
      } else if (roots_size == 1) {
        // One root, if <=, then interval is [r,r]
        if (sgn_condition == LP_SGN_LE_0) {
          lp_interval_construct_point(&x_interval, roots);
          lp_value_destruct(roots);
        } else {
          lp_value_destruct(roots);
          conflict = 1;
          break;
        }
      } else if (roots_size == 2) {
        // Two roots, the interval is either (r0, r1), or [r1, r2]
        int open = sgn_condition == LP_SGN_LT_0;
        lp_interval_construct(&x_interval, roots, open, roots + 1, open);
        lp_value_destruct(roots);
        lp_value_destruct(roots + 1);
      }
      if (trace_is_enabled("polynomial::bounds")) {
        tracef("x_interval = "); lp_interval_print(&x_interval, trace_out); tracef("\n");
      }
      lp_interval_assignment_set_interval(M, x, &x_interval);
      lp_interval_destruct(&x_interval);
    }

    coefficient_destruct(&f);
    coefficient_destruct(&next);
  }

  integer_destruct(&tmp_z);
  integer_destruct(&B_sq);
  integer_destruct(&A4);
  rational_destruct(&D);

  if (ok) {
    if (conflict) {
      return -1;
    } else {
      return 1;
    }
  } else {
    return 0;
  }
}

lp_polynomial_t* lp_polynomial_constraint_explain_infer_bounds(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, int negated, lp_variable_t x) {

  // Same as in infer, just return the polynomial for x

  // Negate the constraint if negated
  if (negated) {
    sgn_condition = lp_sign_condition_negate(sgn_condition);
  }

  const lp_polynomial_context_t* ctx = A->ctx;

  switch (sgn_condition) {
  case LP_SGN_LT_0: // |x| - d < 0 => d < x < d
  case LP_SGN_LE_0: // |x| - d <= 0 => d <= x <= d
    break;
  case LP_SGN_EQ_0: {
    // |x| - d == 0 == both <=, >=
    lp_polynomial_t* p = lp_polynomial_constraint_explain_infer_bounds(A, LP_SGN_LE_0, 0, x);
    if (p) {
      return p;
    }
    return lp_polynomial_constraint_explain_infer_bounds(A, LP_SGN_LE_0, 0, x);
  }
  case LP_SGN_NE_0: // |x| - d != 0 => ?
    return 0;
  case LP_SGN_GT_0: {
    lp_polynomial_t* A_neg = lp_polynomial_new(ctx);
    lp_polynomial_neg(A_neg, A);
    lp_polynomial_t* p = lp_polynomial_constraint_explain_infer_bounds(A_neg, LP_SGN_LT_0, 0, x);
    lp_polynomial_delete(A_neg);
    return p;
  }
  case LP_SGN_GE_0: {
    lp_polynomial_t* A_neg = lp_polynomial_new(ctx);
    lp_polynomial_neg(A_neg, A);
    lp_polynomial_t* p = lp_polynomial_constraint_explain_infer_bounds(A_neg, LP_SGN_LE_0, 0, x);
    lp_polynomial_delete(A_neg);
    return p;
  }
  }

  if (trace_is_enabled("polynomial::bounds")) {
    tracef("lp_polynomial_constraint_explain_infer_bounds("); lp_polynomial_print(A, trace_out); tracef(", "); lp_sign_condition_print(sgn_condition, trace_out); tracef(")\n");
  }

  lp_polynomial_t* result = 0;

  // Make sure we're in the right order
  lp_polynomial_external_clean(A);

  // For each variable in A, see A(x). If
  //
  //    A(x) = Ax^2 + Bx + C,
  //
  // with A and B appropriate constants, then we can write this as
  //
  //   A(x) = a^2 x - 2ab x + b^2 + C - b^2
  //        = (ax - b)^2 + C - b^2
  //
  // If this is the case, and we can recursively rewrite the C further into
  //
  //  (a1 x1 - b1)^2 + ... + (an xn - bn)^2 + C0 - b1^2 - ... - bn^2
  //
  // With D = (b1^2 + ... + bn^2) - C0, and d = sqrt(D), we can infer bounds
  //
  //    (ak x - bk)^2 <= D => -(d + bk)/ak <= x <= (d + bk)/ak
  //
  // In fact, we can just solve
  //
  //    ak^2 - 2akabk x - D^2 == 0
  //
  // and make the interval (r1, r2), [r1], or 0, depending on the number of
  // roots and the sign condition.
  //
  // We can do this easily in the original polynomial, but D is rational so we need
  // to multiply it out.

  lp_integer_t tmp_z, B_sq, A4;
  integer_construct(&tmp_z);
  integer_construct(&B_sq);
  integer_construct(&A4);

  lp_rational_t D;
  rational_construct(&D);

  int ok = 1;
  const coefficient_t* Ak = &A->data;
  while (ok && Ak->type != COEFFICIENT_NUMERIC) {
    // Check if Ak = (ax - b)^2 + Ak-1 = Ax^2 - Bx + ...
    if (trace_is_enabled("polynomial::bounds")) {
      tracef("A_k = "); coefficient_print(ctx, Ak, trace_out); tracef("\n");
      tracef("D = "); rational_print(&D, trace_out); tracef("\n");
    }
    size_t Ak_degree = coefficient_degree(Ak);
    if (Ak_degree == 2) {
      const coefficient_t* A = coefficient_get_coefficient(Ak, 2);
      const coefficient_t* B = coefficient_get_coefficient(Ak, 1);
      if (A->type != COEFFICIENT_NUMERIC || B->type != COEFFICIENT_NUMERIC) {
        ok = 0;
      } else if (integer_sgn(lp_Z, &A->value.num) > 0) {
        //  A(x) = (ax - b)^2 + ...
        //  A = a^2, B=-2ab
        //  A(x) = Ax^2 + Bx + (b^2 + ...)
        //  b = -B/2*a
        //  b^2 = B^2/4*A
        integer_mul(lp_Z, &B_sq, &B->value.num, &B->value.num);
        integer_mul_int(lp_Z, &A4, &A->value.num, 4);
        lp_rational_t tmp_q;
        rational_construct_from_div(&tmp_q, &B_sq, &A4);
        rational_add(&D, &D, &tmp_q);
        rational_destruct(&tmp_q);
        // Go to the next one
        Ak = COEFF(Ak, 0);
      } else {
        // Cannot do square root, must be positive
        ok = 0;
      }
    } else {
      ok = 0;
    }
  }

  if (ok) {
    // D = d^2
    lp_rational_t tmp_q;
    rational_construct_from_integer(&tmp_q, &Ak->value.num);
    rational_sub(&D, &D, &tmp_q);
    rational_destruct(&tmp_q);

    // Construct a polynomial for root finding
    coefficient_t f, next;
    coefficient_construct_copy(ctx, &next, &A->data);
    coefficient_construct(ctx, &f);

    // Now, we traverse again, and get the intervals
    while (next.type != COEFFICIENT_NUMERIC) {
      // Continue to the next one
      coefficient_swap(&f, &next);
      coefficient_assign_int(ctx, &next, 0);
      coefficient_swap(&next, COEFF(&f, 0));
      if (VAR(&f) == x) {
        if (trace_is_enabled("polynomial::bounds")) {
          tracef("f = "); coefficient_print(ctx, &f, trace_out); tracef("\n");
          tracef("D = "); rational_print(&D, trace_out); tracef("\n");
        }
        // Solve Ax^2 + Bx + b^2 - D == 0
        // D is rational p/q we solve
        // B = -2ab => b^2 = B^2/4a^2 = B^2/4A
        // b^2 = B^2/4*A
        const coefficient_t* A = COEFF(&f, 2);
        const coefficient_t* B = COEFF(&f, 1);
        integer_mul(lp_Z, &B_sq, &B->value.num, &B->value.num);
        integer_mul_int(lp_Z, &A4, &A->value.num, 4);
        rational_construct_from_div(&tmp_q, &B_sq, &A4);
        rational_sub(&tmp_q, &tmp_q, &D);
        // Add p/q to polynomial to solve
        const lp_integer_t* p = rational_get_num_ref(&tmp_q);
        const lp_integer_t* q = rational_get_den_ref(&tmp_q);
        coefficient_assign_int(ctx, COEFF(&f, 0), 0);
        coefficient_mul_integer(ctx, &f, &f, q);
        coefficient_assign_integer(ctx, COEFF(&f, 0), p);
        rational_destruct(&tmp_q);
        if (trace_is_enabled("polynomial::bounds")) {
          tracef("f = "); coefficient_print(ctx, &f, trace_out); tracef("\n");
        }
        result = lp_polynomial_new_from_coefficient(ctx, &f);
        break;
      }
    }

    coefficient_destruct(&f);
    coefficient_destruct(&next);
  }

  integer_destruct(&tmp_z);
  integer_destruct(&B_sq);
  integer_destruct(&A4);
  rational_destruct(&D);

  return result;
}

lp_feasibility_set_t* lp_polynomial_root_constraint_get_feasible_set(const lp_polynomial_t* A, size_t root_index, lp_sign_condition_t sgn_condition, int negated, const lp_assignment_t* M) {

  if (trace_is_enabled("polynomial")) {
    tracef("lp_polynomial_root_constraint_get_feasible_set("); lp_polynomial_print(A, trace_out); tracef(", %zu, ", root_index); lp_sign_condition_print(sgn_condition, trace_out); tracef(")\n");
  }

  assert(!lp_polynomial_is_constant(A));

  // Make sure we're in the right order
  lp_polynomial_external_clean(A);

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(A, M, lp_polynomial_top_variable(A));
  }

  // Get the degree of the polynomial, respecting the model
  size_t degree = coefficient_degree_m(A->ctx, &A->data, M);
  if (degree == 0) {
    // Polynomial evaluates to constant, no roots => root_index < 0 is false
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
  assert(roots_size <= degree);

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

int lp_polynomial_check_integrity(const lp_polynomial_t* A) {
  switch (A->data.type) {
  case COEFFICIENT_NUMERIC:
  case COEFFICIENT_POLYNOMIAL:
    return 1;
  default:
    return 0;
  }
}

int lp_polynomial_constraint_resolve_fm(
    const lp_polynomial_t* p1, lp_sign_condition_t p1_sgn,
    const lp_polynomial_t* p2, lp_sign_condition_t p2_sgn,
    const lp_assignment_t* M,
    lp_polynomial_t* R, lp_sign_condition_t* R_sgn,
    lp_polynomial_vector_t* assumptions) {

  lp_polynomial_external_clean(p1);
  lp_polynomial_external_clean(p2);

  lp_variable_t x = lp_polynomial_top_variable(p1);
  if (lp_polynomial_top_variable(p2) != x) {
    return 0;
  }

  if (trace_is_enabled("polynomial::check_input")) {
    check_polynomial_assignment(p1, M, x);
    check_polynomial_assignment(p2, M, x);
  }

  const lp_polynomial_context_t* ctx = p1->ctx;
  assert(p2->ctx == ctx);

  coefficient_t p1_c;
  coefficient_t p2_c;
  coefficient_construct(ctx, &p1_c);
  coefficient_construct(ctx, &p2_c);

  // Reduce, in case we get linear
  coefficient_reductum_m(ctx, &p1_c, &p1->data, M, assumptions);
  coefficient_reductum_m(ctx, &p2_c, &p2->data, M, assumptions);

  int ok = 1;
  if (coefficient_degree(&p1_c) != 1 || coefficient_top_variable(&p1_c) != x) {
    ok = 0;
  }
  if (coefficient_degree(&p2_c) != 1 || coefficient_top_variable(&p2_c) != x) {
    ok = 0;
  }

  if (ok) {

    // Normalize all to be <, <=, ==, or !=
    if (p1_sgn == LP_SGN_GT_0) {
      coefficient_neg(ctx, &p1_c, &p1_c);
      p1_sgn = LP_SGN_LT_0;
    } else if (p1_sgn == LP_SGN_GE_0) {
      coefficient_neg(ctx, &p1_c, &p1_c);
      p1_sgn = LP_SGN_LE_0;
    }
    if (p2_sgn == LP_SGN_GT_0) {
      coefficient_neg(ctx, &p2_c, &p2_c);
      p2_sgn = LP_SGN_LT_0;
    } else if (p2_sgn == LP_SGN_GE_0) {
      coefficient_neg(ctx, &p2_c, &p2_c);
      p2_sgn = LP_SGN_LE_0;
    }

    // Compute the resultant condition
    switch (p1_sgn) {
    case LP_SGN_LT_0:
      switch (p2_sgn) {
      case LP_SGN_LT_0:
        *R_sgn = LP_SGN_LT_0;
        break;
      case LP_SGN_LE_0:
        *R_sgn = LP_SGN_LT_0;
        break;
      case LP_SGN_EQ_0:
        *R_sgn = LP_SGN_LT_0;
        break;
      case LP_SGN_NE_0:
        ok = 0;
        break;
      case LP_SGN_GT_0:
      case LP_SGN_GE_0:
        assert(0);
      }
      break;
    case LP_SGN_LE_0:
      switch (p2_sgn) {
      case LP_SGN_LT_0:
        *R_sgn = LP_SGN_LT_0;
        break;
      case LP_SGN_LE_0:
        *R_sgn = LP_SGN_LE_0;
        break;
      case LP_SGN_EQ_0:
        *R_sgn = LP_SGN_LE_0;
        break;
      case LP_SGN_NE_0:
        ok = 0;
        break;
      case LP_SGN_GT_0:
      case LP_SGN_GE_0:
        assert(0);
      }
      break;
    case LP_SGN_EQ_0:
      switch (p2_sgn) {
      case LP_SGN_LT_0:
      case LP_SGN_LE_0:
      case LP_SGN_EQ_0:
      case LP_SGN_NE_0:
        ok = 0;
        break;
      case LP_SGN_GT_0:
      case LP_SGN_GE_0:
        assert(0);
      }
      break;
    case LP_SGN_NE_0:
      ok = 0;
      break;
    case LP_SGN_GT_0:
    case LP_SGN_GE_0:
      assert(0);
    }

    if (ok) {
      const coefficient_t* p1_lc = coefficient_lc(&p1_c);
      const coefficient_t* p2_lc = coefficient_lc(&p2_c);

      if (p1_lc->type != COEFFICIENT_NUMERIC) {
        lp_polynomial_vector_push_back_coeff(assumptions, p1_lc);
      }
      if (p2_lc->type != COEFFICIENT_NUMERIC) {
        lp_polynomial_vector_push_back_coeff(assumptions, p2_lc);
      }

      int p1_lc_sgn = coefficient_sgn(ctx, p1_lc, M);
      int p2_lc_sgn = coefficient_sgn(ctx, p2_lc, M);

      // The signs must be opposite, unless one of them is ==
      // In that case we can multiply == with negative, still safe
      if (p1_lc_sgn == p2_lc_sgn) {
        if (p1_sgn == LP_SGN_EQ_0) {
          p2_lc_sgn = -p2_lc_sgn;
        } else if (p2_sgn == LP_SGN_EQ_0) {
          p1_lc_sgn = -p1_lc_sgn;
        }
      } else {
        ok = 0;
      }

      if (ok) {
        coefficient_t p1_lc_abs;
        if (p1_lc_sgn > 0) {
          coefficient_construct_copy(ctx, &p1_lc_abs, p1_lc);
        } else {
          coefficient_construct(ctx, &p1_lc_abs);
          coefficient_neg(ctx, &p1_lc_abs, p1_lc);
        }
        coefficient_t p2_lc_abs;
        if (p2_lc_sgn > 0) {
          coefficient_construct_copy(ctx, &p2_lc_abs, p2_lc);
        } else {
          coefficient_construct(ctx, &p2_lc_abs);
          coefficient_neg(ctx, &p2_lc_abs, p2_lc);
        }

//        tracef("p1_c = "); coefficient_print(ctx, &p1_c, trace_out); tracef("\n");
//        tracef("p1_lc_abs = "); coefficient_print(ctx, &p1_lc_abs, trace_out); tracef("\n");
//        tracef("p2_c = "); coefficient_print(ctx, &p2_c, trace_out); tracef("\n");
//        tracef("p2_lc_abs = "); coefficient_print(ctx, &p2_lc_abs, trace_out); tracef("\n");

        coefficient_t R_c;
        coefficient_construct(ctx, &R_c);

        // Compute the resultant polynomial
//        tracef("R = "); coefficient_print(ctx, &R_c, trace_out); tracef("\n");
        coefficient_add_mul(ctx, &R_c, &p1_c, &p2_lc_abs);
//        tracef("R = "); coefficient_print(ctx, &R_c, trace_out); tracef("\n");
        coefficient_add_mul(ctx, &R_c, &p2_c, &p1_lc_abs);
//        tracef("R = "); coefficient_print(ctx, &R_c, trace_out); tracef("\n");

        lp_polynomial_destruct(R);
        lp_polynomial_construct_from_coefficient(R, ctx, &R_c);

        coefficient_destruct(&p1_lc_abs);
        coefficient_destruct(&p2_lc_abs);
        coefficient_destruct(&R_c);
      }
    }
  }

  coefficient_destruct(&p1_c);
  coefficient_destruct(&p2_c);

  return ok;
}

