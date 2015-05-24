/*
 * psc.c
 *
 *  Created on: May 24, 2015
 *      Author: dejan
 */

#include "polynomial/coefficient.h"
#include "polynomial/output.h"

#include "number/integer.h"
#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <assert.h>

STAT_DECLARE(int, coefficient, psc)

/**
 * (non-optimized) Subresultant algorithm, as described in
 *
 * TODO: do the optimized version
 *
 * [2000] Ducos - Optimizations of the subresultant algorithm.
 */
void coefficient_psc(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* P, const coefficient_t* Q) {

  TRACE("coefficient", "coefficient_psc()\n");
  STAT(coefficient, psc) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
  }

  assert(P->type == COEFFICIENT_POLYNOMIAL);
  assert(Q->type == COEFFICIENT_POLYNOMIAL);

  lp_variable_t x = VAR(P);
  assert(VAR(Q) == x);

  size_t P_deg = coefficient_degree(P);
  size_t Q_deg = coefficient_degree(Q);
  assert(P_deg >= Q_deg);

  // S = []
  int S_size = 0;

  // s = lc(Q)^(deg(P) - deg(Q)
  coefficient_t s;
  coefficient_construct(ctx, &s);
  coefficient_pow(ctx, &s, coefficient_lc(Q), P_deg - Q_deg);

  // Set the final position
  coefficient_assign(ctx, S + Q_deg, &s);

  // A = Q, B = prem(P, -Q)
  coefficient_t A, B;
  coefficient_construct_copy(ctx, &A, Q);
  coefficient_construct_copy(ctx, &B, Q);
  coefficient_neg(ctx, &B, &B);
  coefficient_prem(ctx, &B, P, &B);

  // Some temporaries
  coefficient_t C, pow;
  coefficient_construct(ctx, &C);
  coefficient_construct(ctx, &pow);

  for (;;) {

    if (trace_is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
    }

    // d = deg(A); e = deg(B)
    size_t d = coefficient_degree_safe(ctx, &A, x);
    size_t e = coefficient_degree_safe(ctx, &B, x);
    assert(d > e);

    // Holds:
    //   A ~ S_d   if d = deg(Q)
    //   A = S_d   if d < deg(Q)
    //   B = S_d-1, s = lc(S_d) for d <= deg(Q)

    if (coefficient_is_zero(ctx, &B)) {
      break;
    }

    // S = [B; S]
    if (coefficient_degree_safe(ctx, &B, x) == d - 1) {
      coefficient_assign(ctx, S + S_size, coefficient_lc_safe(ctx, &B, x));
    }
    S_size ++;
    if (trace_is_enabled("coefficient::resultant")) {
      tracef("S[%d] = ", S_size - 1); coefficient_print(ctx, S + S_size - 1, trace_out); tracef("\n");
    }

    // Holds:
    //   S = [S_d-1, S_d, ...]

    int delta = d - e;
    if (delta > 1) {
      // C = (lc(B)^(delta-1)*B)/(s^(delta-1))
      coefficient_pow(ctx, &pow, coefficient_lc_safe(ctx, &B, x), delta-1);
      coefficient_mul(ctx, &C, &pow, &B);
      coefficient_pow(ctx, &pow, &s, delta-1);
      coefficient_div(ctx, &C, &C, &pow);
      // S = [C; S]
      if (coefficient_degree_safe(ctx, &C, x) == e) {
        coefficient_assign(ctx, S + S_size, coefficient_lc_safe(ctx, &C, x));
      }
      S_size ++;
      if (trace_is_enabled("coefficient::resultant")) {
        tracef("S[%d] = ", S_size - 1); coefficient_print(ctx, S + S_size - 1, trace_out); tracef("\n");
      }
    } else {
      // C = B
      coefficient_assign(ctx, &C, &B);
    }

    // Holds:
    //   C = S_e, S = [S_e, ...]

    if (e == 0) {
      break;
    }

    if (trace_is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("lc(A) = "); coefficient_print(ctx, coefficient_lc_safe(ctx, &A, x), trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
      tracef("s = "); coefficient_print(ctx, &s, trace_out); tracef("\n");
      tracef("delta = %d\n", delta);
    }

    // B = prem(A, -B)/(s^delta*lc(A))
    coefficient_neg(ctx, &B, &B);
    coefficient_prem(ctx, &B, &A, &B);
    coefficient_pow(ctx, &pow, &s, delta);
    coefficient_mul(ctx, &pow, &pow, coefficient_lc_safe(ctx, &A, x));
    coefficient_div(ctx, &B, &B, &pow);

    if (trace_is_enabled("coefficient::resultant")) {
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
    }

    // Holds:
    //   B = S_e-1

    if (trace_is_enabled("coefficient::resultant")) {
      tracef("C = "); coefficient_print(ctx, &C, trace_out); tracef("\n");
      tracef("lc(A) = "); coefficient_print(ctx, coefficient_lc_safe(ctx, &A, x), trace_out); tracef("\n");
    }

    // A = C, s = lc(A)
    coefficient_swap(&A, &C);
    coefficient_assign(ctx, &s, coefficient_lc_safe(ctx, &A, x));

    if (trace_is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("s = "); coefficient_print(ctx, &s, trace_out); tracef("\n");
    }
  }

  // Reverse S
  int i = 0, j = S_size - 1;
  while (i < j) {
    coefficient_swap(S + i, S + j);
    i ++;
    j --;
  }

  // Remove temps
  coefficient_destruct(&A);
  coefficient_destruct(&B);
  coefficient_destruct(&C);
  coefficient_destruct(&pow);
  coefficient_destruct(&s);
}

