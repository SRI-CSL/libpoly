/*
 * polynomial_internal.c
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#include "polynomial/internal.h"

#include "utils/debug_trace.h"

#include <assert.h>
#include <malloc.h>

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

void polynomial_external_clean(const polynomial_t* A_const) {
  if (A_const->flags.external && !coefficient_ops.in_order(A_const->ctx, &A_const->data)) {
    polynomial_t* A = (polynomial_t*) A_const;
    coefficient_ops.order(A->ctx, &A->data);
    A->flags.primitive = 0;
  }
}

inline
void polynomial_set_context(polynomial_t* A, const polynomial_context_t* ctx) {
  if (A->ctx != ctx) {
    if (A->ctx && A->flags.external) {
      polynomial_context_ops.detach((polynomial_context_t*)A->ctx);
    }
    A->ctx = ctx;
    if (A->ctx && A->flags.external) {
      polynomial_context_ops.attach((polynomial_context_t*)A->ctx);
    }
  }
}

void polynomial_construct(polynomial_t* A, const polynomial_context_t* ctx, int external) {
  A->ctx = 0;
  A->flags.external = external;
  A->flags.prime = 1;
  A->flags.primitive = 1;
  A->flags.univariate = 1;
  polynomial_set_context(A, ctx);
  coefficient_ops.construct(ctx, &A->data);
}

void polynomial_construct_from_coefficient(polynomial_t* A, const polynomial_context_t* ctx,
    const coefficient_t* from, int external) {
  A->ctx = 0;
  A->flags.prime = 0;
  A->flags.primitive = 0;
  A->flags.univariate = 0;
  A->flags.external = external;
  polynomial_set_context(A, ctx);
  coefficient_ops.construct_copy(A->ctx, &A->data, from);
}

void polynomial_construct_copy(polynomial_t* A, const polynomial_t* from, int external) {
  A->ctx = 0;
  A->flags = from->flags;
  A->flags.external = external;
  polynomial_set_context(A, from->ctx);
  coefficient_ops.construct_copy(A->ctx, &A->data, &from->data);
}

/** Construct a simple polynomial c*x^n */
void polynomial_construct_simple(
    polynomial_t* A, const polynomial_context_t* ctx, int external,
    const integer_t* c, variable_t x, unsigned n)
{
  A->ctx = 0;
  A->flags.external = external;
  A->flags.prime = 0;
  A->flags.primitive = 0;
  A->flags.univariate = 1;
  polynomial_set_context(A, ctx);
  coefficient_ops.construct_simple(ctx, &A->data, c, x, n);
}

void polynomial_destruct(polynomial_t* A) {
  coefficient_ops.destruct(&A->data);
  if (A->flags.external) {
    polynomial_context_ops.detach((polynomial_context_t*)A->ctx);
  }
}

polynomial_t* polynomial_alloc(void) {
  polynomial_t* new = malloc(sizeof(polynomial_t));
  return new;
}

polynomial_t* polynomial_new(const polynomial_context_t* ctx, int external) {
  polynomial_t* new = polynomial_alloc();
  polynomial_construct(new, ctx, external);
  return new;
}

#define SWAP(type, x, y) { type tmp = x; x = y; y = tmp; }

void polynomial_swap(polynomial_t* A1, polynomial_t* A2) {
  // Swap everything, but keep the external flags
  polynomial_t tmp = *A1; *A1 = *A2; *A2 = tmp;
  SWAP(unsigned, A1->flags.external, A2->flags.external);
}

void polynomial_assign(polynomial_t* A, const polynomial_t* from) {
  if (A != from) {
    polynomial_set_context(A, from->ctx);
    if (A->flags.external) {
      A->flags = from->flags;
      A->flags.external = 1;
    } else {
      A->flags = from->flags;
      A->flags.external = 0;
    }
    coefficient_ops.assign(A->ctx, &A->data, &from->data);
  }
}

const polynomial_context_t* polynomial_context(const polynomial_t* A) {
  return A->ctx;
}

variable_t polynomial_top_variable(const polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_ops.top_variable(&A->data);
}

size_t polynomial_degree(const polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_ops.degree(&A->data);
}

void polynomial_get_coefficient(polynomial_t* C_p, const polynomial_t* A, size_t k) {
  polynomial_external_clean(A);

  if (k > polynomial_degree(A)) {
    polynomial_t result;
    polynomial_construct(&result, A->ctx, 0);
    polynomial_swap(C_p, &result);
    polynomial_destruct(&result);
  } else {
    const coefficient_t* C = coefficient_ops.get_coefficient(&A->data, k);
    polynomial_t result;
    polynomial_construct_from_coefficient(&result, A->ctx, C, 0);
    polynomial_swap(C_p, &result);
    polynomial_destruct(&result);
  }
}

int polynomial_is_constant(const polynomial_t* A) {
  return coefficient_ops.is_constant(&A->data);
}

int polynomial_is_zero(const polynomial_t* A) {
  return coefficient_ops.is_zero(A->ctx, &A->data);
}

int polynomial_sgn(const polynomial_t* A, const assignment_t* m) {
  return coefficient_ops.sgn(A->ctx, &A->data, m);
}

int polynomial_cmp(const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_cmp(%P, %P)\n", A1, A2);

  if (!polynomial_context_ops.equal(A1->ctx, A2->ctx)) {
    // random order for different contexts
    return A1 - A2;
  }
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  int cmp = coefficient_ops.cmp(A1->ctx, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_cmp(%P, %P) => %d\n", A1, A2, cmp);

  return cmp;
}

int polynomial_cmp_type(const polynomial_t* A1, const polynomial_t* A2) {
  const polynomial_context_t* ctx = A1->ctx;
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  return coefficient_ops.cmp_type(ctx, &A1->data, &A2->data);
}

int polynomial_divides(const polynomial_t* A1, const polynomial_t* A2) {
  if (!polynomial_context_ops.equal(A1->ctx, A2->ctx)) {
    return 0;
  }
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);
  return coefficient_ops.divides(A1->ctx, &A1->data, &A2->data);
}

int polynomial_print(const polynomial_t* A, FILE* out) {
  polynomial_external_clean(A);
  return coefficient_ops.print(A->ctx, &A->data, out);
}

char* polynomial_to_string(const polynomial_t* A) {
  polynomial_external_clean(A);
  return coefficient_ops.to_string(A->ctx, &A->data);
}

inline
int polynomial_same_var_univariate(const polynomial_t* A1, const polynomial_t* A2) {
  if (!A1->flags.univariate) {
    return 0;
  }
  if (!A2->flags.univariate) {
    return 0;
  }
  if (A1->data.type == COEFFICIENT_NUMERIC) {
    return 1;
  }
  if (A2->data.type == COEFFICIENT_NUMERIC) {
    return 1;
  }
  // Both univariate and not constant
  return A1->data.value.rec.x == A2->data.value.rec.x;
}

void polynomial_add(polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_add(%P, %P, %P)\n", S, A1, A2);

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  S->flags.prime = 0;
  S->flags.primitive = 0;
  S->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(S, A1->ctx);

  coefficient_ops.add(S->ctx, &S->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_add() => %P\n", S);
}

void polynomial_sub(polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2) {
  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  TRACE("polynomial", "polynomial_sub(%P, %P, %P)\n", S, A1, A2);

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  S->flags.prime = 0;
  S->flags.primitive = 0;
  S->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(S, A1->ctx);

  coefficient_ops.sub(S->ctx, &S->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_sub() => %P\n", S);
}

void polynomial_neg(polynomial_t* N, const polynomial_t* A) {

  polynomial_external_clean(A);

  N->flags.prime = 0;
  N->flags.primitive = 0;
  N->flags.univariate = A->flags.univariate;

  polynomial_set_context(N, N->ctx);

  coefficient_ops.neg(N->ctx, &N->data, &A->data);
}

void polynomial_mul(polynomial_t* P, const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_mul(%P, %P, %P)\n", P, A1, A2);

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  P->flags.prime = 0;
  P->flags.primitive = 0;
  P->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(P, A1->ctx);

  coefficient_ops.mul(P->ctx, &P->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_mul() => %P\n", P);
}

void polynomial_shl(polynomial_t* S, const polynomial_t* A, unsigned n) {

  polynomial_external_clean(A);

  S->flags.prime = 0;
  S->flags.primitive = 0;
  S->flags.univariate = A->flags.univariate;

  polynomial_set_context(S, A->ctx);

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  coefficient_ops.shl(S->ctx, &S->data, &A->data, VAR(&A->data), n);
}

void polynomial_pow(polynomial_t* P, const polynomial_t* A, unsigned n) {

  TRACE("polynomial", "polynomial_pow(%P, %P)\n", P, A);

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);

  P->flags.prime = 0;
  P->flags.primitive = 0;
  P->flags.univariate = A->flags.univariate;

  polynomial_set_context(P, A->ctx);

  coefficient_ops.pow(P->ctx, &P->data, &A->data, n);

  TRACE("polynomial", "polynomial_pow() => %P\n", P);
}

void polynomial_add_mul(polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2) {
  const polynomial_context_t* ctx = A1->ctx;

  assert(polynomial_context_ops.equal(S->ctx, ctx));
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));

  polynomial_external_clean(S);
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  S->flags.prime = 0;
  S->flags.primitive = 0;
  S->flags.univariate = polynomial_same_var_univariate(A1, A2) && polynomial_same_var_univariate(S, A1);

  coefficient_ops.add_mul(ctx, &S->data, &A1->data, &A2->data);
}

void polynomial_sub_mul(polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2) {
  const polynomial_context_t* ctx = A1->ctx;

  assert(polynomial_context_ops.equal(S->ctx, ctx));
  assert(polynomial_context_ops.equal(A1->ctx, ctx));
  assert(polynomial_context_ops.equal(A2->ctx, ctx));

  polynomial_external_clean(S);
  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  S->flags.prime = 0;
  S->flags.primitive = 0;
  S->flags.univariate = polynomial_same_var_univariate(A1, A2) && polynomial_same_var_univariate(S, A1);

  coefficient_ops.sub_mul(ctx, &S->data, &A1->data, &A2->data);
}

void polynomial_div(polynomial_t* D, const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_div(%P, %P, %P)\n", D, A1, A2);

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  D->flags.prime = 0;
  D->flags.primitive = 0;
  D->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(D, A1->ctx);

  coefficient_ops.div(D->ctx, &D->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_div() => %P\n", D);
}

void polynomial_rem(polynomial_t* R, const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_rem(%P, %P, %P)\n", R, A1, A2);

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  R->flags.prime = 0;
  R->flags.primitive = 0;
  R->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(R, A1->ctx);

  coefficient_ops.rem(R->ctx, &R->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_rem() => %P\n", R);
}

void polynomial_divrem(polynomial_t* D, polynomial_t* R, const polynomial_t* A1, const polynomial_t* A2) {

  TRACE("polynomial", "polynomial_divrem(%P, %P, %P, %P)\n", D, R, A1, A2);

  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  R->flags.prime = 0;
  R->flags.primitive = 0;
  R->flags.univariate = polynomial_same_var_univariate(A1, A2);

  polynomial_set_context(D, A1->ctx);
  polynomial_set_context(R, A1->ctx);

  coefficient_ops.divrem(D->ctx, &R->data, &R->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_rem() => (%P, %P)\n", D, R);
}

void polynomial_derivative(polynomial_t* A_d, const polynomial_t* A) {

  TRACE("polynomial", "polynomial_derivative(%P, %P)\n", A_d, A);

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);

  // primitive polynomials are not invariant with derivative
  // A = x*2 + 2x
  // A' = 2x + 2

  A_d->flags.prime = 0;
  A_d->flags.primitive = 0;
  A_d->flags.univariate = A->flags.univariate;

  polynomial_set_context(A_d, A->ctx);

  coefficient_ops.derivative(A_d->ctx, &A_d->data, &A->data);

  TRACE("polynomial", "polynomial_derivative() => %P\n", A_d);
}

void polynomial_gcd(polynomial_t* gcd, const polynomial_t* A1, const polynomial_t* A2) {
  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  TRACE("polynomial", "polynomial_gcd(%P, %P)\n", A1, A2);

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A1->ctx->var_order, A1->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(gcd, A1->ctx);

  gcd->flags.prime = 0;
  gcd->flags.primitive = A1->flags.primitive && A2->flags.primitive;
  gcd->flags.univariate = polynomial_same_var_univariate(A1, A2);

  coefficient_ops.gcd(gcd->ctx, &gcd->data, &A1->data, &A2->data);

  TRACE("polynomial", "polynomial_gcd() => %P\n", gcd);
}

void polynomial_lcm(polynomial_t* lcm, const polynomial_t* A1, const polynomial_t* A2) {
  assert(polynomial_context_ops.equal(A1->ctx, A2->ctx));

  polynomial_external_clean(A1);
  polynomial_external_clean(A2);

  polynomial_set_context(lcm, A1->ctx);

  lcm->flags.prime = 0;
  lcm->flags.primitive = A1->flags.primitive && A2->flags.primitive;
  lcm->flags.univariate = polynomial_same_var_univariate(A1, A2);

  coefficient_ops.lcm(lcm->ctx, &lcm->data, &A1->data, &A2->data);
}

void polynomial_reduce(
    const polynomial_t* A, const polynomial_t* B,
    polynomial_t* P, polynomial_t* Q, polynomial_t* R) {

  TRACE("polynomial", "polynomial_reduce(%P, %P)\n", A, B);

  const polynomial_context_t* ctx = A->ctx;

  assert(polynomial_context_ops.equal(B->ctx, ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  polynomial_set_context(P, ctx);
  polynomial_set_context(Q, ctx);
  polynomial_set_context(R, ctx);

  coefficient_ops.reduce(ctx, &A->data, &B->data, &P->data, &Q->data, &R->data, 1);

  TRACE("polynomial", "polynomial_derivative() =>\n");
  TRACE("polynomial", "\t P = %P\n", P);
  TRACE("polynomial", "\t Q = %P\n", Q);
  TRACE("polynomial", "\t R = %P\n", R);
}

void polynomial_psc(polynomial_t** psc, const polynomial_t* A, const polynomial_t* B) {

  TRACE("polynomial", "polynomial_psc(%P, %P)\n", A, B);

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  size_t A_deg = polynomial_degree(A);
  size_t B_deg = polynomial_degree(B);

  if (A_deg < B_deg) {
    polynomial_psc(psc, B, A);
    return;
  }

  const polynomial_context_t* ctx = A->ctx;
  assert(polynomial_context_ops.equal(B->ctx, ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  // Allocate the space for the result
  size_t size = B_deg + 1;
  coefficient_t* psc_coeff = malloc(sizeof(coefficient_t)*size);
  int i;
  for (i = 0; i < size; ++ i) {
    coefficient_ops.construct(ctx, psc_coeff + i);
  }

  // Compute
  coefficient_ops.psc(ctx, psc_coeff, &A->data, &B->data);

  // Construct the output (one less, we ignore the final 1)
  for (i = 0; i < size; ++ i) {
    polynomial_t tmp;
    polynomial_construct_from_coefficient(&tmp, ctx, psc_coeff + i, 0);
    polynomial_swap(&tmp, psc[i]);
    polynomial_destruct(&tmp);
    coefficient_ops.destruct(&psc_coeff[i]);
  }

  free(psc_coeff);

  if (debug_trace_ops.is_enabled("polynomial")) {
    for (i = 0; i < size; ++ i) {
      tracef("PSC[%d] = %P\n", i, psc[i]);
    }
  }
}

void polynomial_resultant(polynomial_t* res, const polynomial_t* A, const polynomial_t* B) {

  TRACE("polynomial", "polynomial_resultant(%P, %P)\n", A, B);

  assert(A->data.type == COEFFICIENT_POLYNOMIAL);
  assert(B->data.type == COEFFICIENT_POLYNOMIAL);
  assert(VAR(&A->data) == VAR(&B->data));

  const polynomial_context_t* ctx = A->ctx;
  assert(polynomial_context_ops.equal(B->ctx, ctx));

  if (debug_trace_ops.is_enabled("polynomial")) {
    variable_order_simple_ops.print((variable_order_simple_t*) A->ctx->var_order, A->ctx->var_db, trace_out);
    tracef("\n");
  }

  polynomial_external_clean(A);
  polynomial_external_clean(B);

  // Compute
  coefficient_ops.resultant(ctx, &res->data, &A->data, &B->data);

  if (debug_trace_ops.is_enabled("polynomial")) {
    tracef("polynomial_resultant(%P, %P) => %P\n", A, B, res);
  }
}

