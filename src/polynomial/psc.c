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
 * [2000] Ducos - Optimizations of the subresultant algorithm.
 */
void coefficient_psc_unoptimized(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* P, const coefficient_t* Q) {

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

static
void S_e_optimized(const lp_polynomial_context_t* ctx, const coefficient_t* S_d, const coefficient_t* S_d_1, coefficient_t* S_e, lp_variable_t top_var) {
  // n = deg(S_d) - deg(S_d_1) - 1
  size_t n = coefficient_degree_safe(ctx, S_d, top_var) - coefficient_degree_safe(ctx, S_d_1, top_var) - 1;
  // if n == 0 then return S_d_1
  if (n == 0) {
    coefficient_assign(ctx, S_e, S_d_1);
    return;
  }
  // (x, y) = (lc(S_d-1), lc(S_d))
  coefficient_t x, y;
  coefficient_construct_copy(ctx, &x, S_d_1);
  coefficient_construct_copy(ctx, &y, S_d);
  // a = 2^floor(log2(n)), a <= n < 2*a
  size_t a = 1;
  while ((a >> 1) < n) { a >>= 1; }
  // c = x
  coefficient_t c;
  coefficient_construct_copy(ctx, &c, &x);
  // n = n - a
  n = n - a;
  // loop
  for (;;) {
    // exit when a = 1
    if (a == 1) { break; }
    // a = a/2, c = c^2/y
    a = a / 2;
    coefficient_mul(ctx, &c, &c, &c);
    coefficient_div(ctx, &c, &c, &y);
    // if n >= a then c = (c*x)/y; n = n - a
    if (n >= a) {
      coefficient_mul(ctx, &c, &c, &x);
      coefficient_div(ctx, &c, &c, &y);
      n = n - a;
    }
  }
  // Return (c*S_d-1)/y
  coefficient_mul(ctx, &c, &c, S_d_1);
  coefficient_div(ctx, S_e, &c, &y);
  // Remove temps
  coefficient_destruct(&c);
  coefficient_destruct(&x);
  coefficient_destruct(&y);
}

static
void S_e_1_optimized(const lp_polynomial_context_t* ctx, const coefficient_t* A, const coefficient_t* S_d_1, const coefficient_t* S_e, const coefficient_t* s_d, coefficient_t* S_e_1, lp_variable_t top_var) {
  // (d, e) = (deg(A), deg(S_d-1))
  size_t d = coefficient_degree_safe(ctx, A, top_var);
  size_t e = coefficient_degree_safe(ctx, S_d_1, top_var);
  // (c_d-1, s_e) = (lc(S_d-1), lc(S_e)
  const coefficient_t* c_d_1;
  const coefficient_t* s_e;
  c_d_1 = coefficient_lc_safe(ctx, S_d_1, top_var);
  s_e = coefficient_lc_safe(ctx, S_e, top_var);
  // for j in 0...e-1 loop
  //   H_j = s_e X^j
  assert(e > 0);
  size_t j;
  coefficient_t* H = malloc(sizeof(coefficient_t)*d);
  assert(e < d);
  for (j = 0; j < e; j ++) {
    coefficient_construct_copy(ctx, H + j, s_e);
    coefficient_shl(ctx, H + j, H + j, top_var, j);
  }
  // H_e  = s_e * X^e - S_e
  coefficient_construct_copy(ctx, H + e, s_e);
  coefficient_shl(ctx, H + e, H + e, top_var, e);
  coefficient_sub(ctx, H + e, H + e, S_e);
  // for j in e + 1 ... d - 1
  //    H_j = x*H_j-1 - pi_e(X*H_j-1)*S_d-1/c_d-1
  coefficient_t tmp;
  coefficient_construct(ctx, &tmp);
  for (j = e + 1; j < d; ++ j) {
    // H_j = x*H_j-1
    coefficient_construct_copy(ctx, H + j, H + j - 1);
    coefficient_shl(ctx, H + j, H + j, top_var, 1);
    // pi_e = e-th coefficient
    const coefficient_t* pi_e = coefficient_get_coefficient_safe(ctx, H + j, e, top_var);
    // second term
    coefficient_mul(ctx, &tmp, pi_e, S_d_1);
    coefficient_div(ctx, &tmp, &tmp, c_d_1);
    // Do the subtraction
    coefficient_sub(ctx, H + j, H + j, &tmp);
  }
  // D = sum{j < d} pi_j(A)*H_j / lc(A)
  coefficient_t D;
  coefficient_construct(ctx, &D);
  for (j = 0; j < d; ++ j) {
    coefficient_add_mul(ctx, &D, coefficient_get_coefficient_safe(ctx, A, j, top_var), H + j);
  }
  coefficient_div(ctx, &D, &D, coefficient_lc_safe(ctx, A, top_var));
  // Finial result = (-1)^(d-e+1) * (c_d-1)(x*H_d-1+D)-pi_e(x*H_d-1)S_d-1)/s_d
  coefficient_t result;
  coefficient_construct(ctx, &result);
  coefficient_shl(ctx, &tmp, H + d - 1, top_var, 1);
  coefficient_mul(ctx, &result, coefficient_get_coefficient_safe(ctx, &tmp, e, top_var), S_d_1);
  coefficient_add(ctx, &tmp, &tmp, &D);
  coefficient_mul(ctx, &tmp, &tmp, c_d_1);
  coefficient_sub(ctx, &result, &tmp, &result);
  coefficient_div(ctx, &result, &result, s_d);
  if ((d - e + 1) % 2) {
    coefficient_neg(ctx, &result, &result);
  }
  coefficient_swap(&result, S_e_1);
  // Remove temps
  coefficient_destruct(&result);
  coefficient_destruct(&D);
  coefficient_destruct(&tmp);
  for (j = 0; j < d; ++ j) {
    coefficient_destruct(H + j);
  }
  free(H);
}


/**
 * (optimized) Subresultant algorithm, as described in
 *
 * [2000] Ducos - Optimizations of the subresultant algorithm.
 */
void coefficient_psc_optimized(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* P, const coefficient_t* Q) {

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
      // Optimized calculation of S_e into C
      S_e_optimized(ctx, S+d, S+d-1, &C, x);
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

    // B = Optimized calculation of S_e-1
    // at this point C = S_e, B = S_d-1,
    S_e_1_optimized(ctx, &A, &B, &C, &s, &B, x);

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

void coefficient_psc(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* P, const coefficient_t* Q) {
  coefficient_psc_unoptimized(ctx, S, P, Q);
}

