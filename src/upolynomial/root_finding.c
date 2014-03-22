/*
 * univariate_polynomial_root_finding.c
 *
 *  Created on: Jan 13, 2014
 *      Author: dejan
 */

#include "upolynomial/root_finding.h"
#include "upolynomial/upolynomial_dense.h"

#include <assert.h>
#include <malloc.h>
#include <stdlib.h>

#include "utils/debug_trace.h"

/** Negative infinity */
#define INF_N (void*) 0
/** Positive infinity */
#define INF_P (void*) 1

void upolynomial_compute_sturm_sequence(const upolynomial_t* f, upolynomial_dense_t* S, size_t* size) {

  TRACE("roots", "upolynomial_compute_sturm_sequence(%P)\n", f);

  integer_t a;
  integer_ops.construct_from_int(Z, &a, 0);

  // Min size for the polynomials
  size_t f_deg = upolynomial_ops.degree(f);

  // f[0] = pp(f)
  upolynomial_dense_ops.construct_p(&S[0], f_deg + 1, f);
  upolynomial_dense_ops.mk_primitive_Z(&S[0], 1);

  if (debug_trace_ops.is_enabled("roots")) {
    tracef("S[0] = ");
    upolynomial_dense_ops.print(&S[0], trace_out);
    tracef("\n");
  }

  // f[1] = f'
  upolynomial_dense_ops.construct(&S[1], f_deg + 1);
  upolynomial_dense_ops.derivative(Z, &S[0], &S[1]);
  upolynomial_dense_ops.mk_primitive_Z(&S[1], 1);

  if (debug_trace_ops.is_enabled("roots")) {
    tracef("S[1] = ");
    upolynomial_dense_ops.print(&S[1], trace_out);
    tracef("\n");
  }

  // Until we hit a constant polynomial
  int i = 1;
  while (S[i].size > 1) {
    i ++;
    upolynomial_dense_ops.construct(&S[i], f_deg + 1);
    // Compute a*S[i-2] = div*S[i-1] + b*S[i]
    upolynomial_dense_ops.reduce_Z(&S[i-2], &S[i-1], &a, &S[i]);
    upolynomial_dense_ops.mk_primitive_Z(&S[i], 0);

    // If the coefficient of the reduction is not negative, we have to flip the
    // sign to get a Sturm sequence
    if (integer_ops.sgn(Z, &a) > 0) {
      upolynomial_dense_ops.negate(&S[i], Z);
    }

    if (debug_trace_ops.is_enabled("roots")) {
      tracef("S[%d] = ", i); upolynomial_dense_ops.print(&S[i], trace_out); tracef("\n");
    }
  }

  // Remove the temp
  integer_ops.destruct(&a);

  // Size
  *size = i + 1;
}

/**
 * Compute the number sgn_changes(a) given Sturm sequence and a. If a is 0 or
 * 1 as a pointer, we evaluate at -inf, +inf respectively.
 */
int sturm_seqence_count_sign_changes(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const rational_t* a, int max_changes)
{
  int i, sgn_a, sgn_a_previous = 0, sgn_a_changes_count = 0;
  for (i = 0; i < sturm_sequence_size && sgn_a_changes_count < max_changes; ++ i) {
    // Get the sign of S[i] at a
    sgn_a = (a == INF_N) ? upolynomial_dense_ops.sgn_at_minus_inf(&sturm_sequence[i])
          : (a == INF_P) ? upolynomial_dense_ops.sgn_at_plus_inf(&sturm_sequence[i])
          : upolynomial_dense_ops.sgn_at_rational(&sturm_sequence[i], a);
    // Compare with the previous non-zero sign
    if (sgn_a_previous == 0) {
      sgn_a_previous = sgn_a;
    } else if (sgn_a && sgn_a * sgn_a_previous < 0) {
      sgn_a_changes_count++;
      sgn_a_previous = sgn_a;
    }
  }
  return sgn_a_changes_count;
}


/**
 * Compute the number sgn_changes(a) given Sturm sequence and a. If a is 0 or
 * 1 as a pointer, we evaluate at -inf, +inf respectively.
 */
int sturm_seqence_count_sign_changes_dyadic(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const dyadic_rational_t* a, int max_changes)
{
  int i, sgn_a, sgn_a_previous = 0, sgn_a_changes_count = 0;
  for (i = 0; i < sturm_sequence_size && sgn_a_changes_count < max_changes; ++ i) {
    // Get the sign of S[i] at a
    sgn_a = (a == INF_N) ? upolynomial_dense_ops.sgn_at_minus_inf(&sturm_sequence[i])
          : (a == INF_P) ? upolynomial_dense_ops.sgn_at_plus_inf(&sturm_sequence[i])
          : upolynomial_dense_ops.sgn_at_dyadic_rational(&sturm_sequence[i], a);
    // Compare with the previous non-zero sign
    if (sgn_a_previous == 0) {
      sgn_a_previous = sgn_a;
    } else if (sgn_a && sgn_a * sgn_a_previous < 0) {
      sgn_a_changes_count++;
      sgn_a_previous = sgn_a;
    }
  }
  return sgn_a_changes_count;
}

/** Count roots in the given interval */
int sturm_seqence_count_roots(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const interval_t* interval)
{
  // Sign changes at a (or -inf)
  int a_sgn_changes =
      interval ? sturm_seqence_count_sign_changes(sturm_sequence, sturm_sequence_size, &interval->a, sturm_sequence_size)
               : sturm_seqence_count_sign_changes(sturm_sequence, sturm_sequence_size, INF_N, sturm_sequence_size);
  // Sign changes a b (or +inf)
  int b_sgn_changes =
      interval ? sturm_seqence_count_sign_changes(sturm_sequence, sturm_sequence_size, &interval->b, sturm_sequence_size)
               : sturm_seqence_count_sign_changes(sturm_sequence, sturm_sequence_size, INF_P, sturm_sequence_size);
  // Number of roots in (a, b]
  int root_count = a_sgn_changes - b_sgn_changes;
  // Adjust for the possible zeros at a dn b
  if (interval) {
    if (interval->b_open && upolynomial_dense_ops.sgn_at_rational(&sturm_sequence[0], &interval->b) == 0) {
      root_count --;
    }
    if (!interval->a_open && upolynomial_dense_ops.sgn_at_rational(&sturm_sequence[0], &interval->a) != 0) {
      root_count ++;
    }
  }
  // Done
  return root_count;
}

/** Count roots in the given interval */
int sturm_seqence_count_roots_dyadic(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const dyadic_interval_t* interval)
{
  // Sign changes at a (or -inf)
  int a_sgn_changes =
      interval ? sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval->a, sturm_sequence_size)
               : sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, INF_N, sturm_sequence_size);
  // Sign changes a b (or +inf)
  int b_sgn_changes =
      interval ? sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval->b, sturm_sequence_size)
               : sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, INF_P, sturm_sequence_size);
  // Number of roots in (a, b]
  int root_count = a_sgn_changes - b_sgn_changes;
  // Adjust for the possible zeros at a dn b
  if (interval) {
    if (interval->b_open && upolynomial_dense_ops.sgn_at_dyadic_rational(&sturm_sequence[0], &interval->b) == 0) {
      root_count --;
    }
    if (!interval->a_open && upolynomial_dense_ops.sgn_at_dyadic_rational(&sturm_sequence[0], &interval->a) != 0) {
      root_count ++;
    }
  }
  // Done
  return root_count;
}

/**
 * Count the distinct real roots using the Sturm's sequence. For a square-free
 * polynomial p, we can compute the number of distinct real roots in the
 * interval (a, b]. We therefore do a square-free decomposition, check if b is
 * a zero and then count the roots in (a, b].
 */
int upolynomial_roots_count_sturm(const upolynomial_t* f, const interval_t* interval) {

  assert(f->K == Z);

  if (debug_trace_ops.is_enabled("roots")) {
    tracef("upolynomial_root_count_sturm(%P, ", f); interval_ops.print(interval, trace_out); tracef("\n");
  }

  int total_count = 0;

  // Special case for the constants
  if (upolynomial_ops.degree(f) == 0) {
    assert(!upolynomial_ops.is_zero(f));
    return 0;
  }

  // Get the square-free factorization and then count roots for each factor.
  upolynomial_factors_t* square_free_factors = upolynomial_ops.factor_square_free(f);

  int factor_i;
  for (factor_i = 0; factor_i < square_free_factors->size; ++ factor_i) {
    // The factor we are working with
    const upolynomial_t* factor = square_free_factors->factors[factor_i];

    // Compute the Sturm sequence for the factor
    size_t sturm_sequence_size;
    upolynomial_dense_t* sturm_sequence = malloc((upolynomial_ops.degree(factor) + 1)*sizeof(upolynomial_dense_t));
    upolynomial_compute_sturm_sequence(factor, sturm_sequence, &sturm_sequence_size);

    // Add the number of roots
    total_count += sturm_seqence_count_roots(sturm_sequence, sturm_sequence_size, interval);

    // Destroy the temporaries
    int i;
    for (i = 0; i < sturm_sequence_size; ++i) {
      upolynomial_dense_ops.destruct(sturm_sequence + i);
    }
    free(sturm_sequence);
  }

  // Destroy the factors
  upolynomial_factors_ops.destruct(square_free_factors, 1);

  // Return the total number of roots
  return total_count;
}

/**
 * Recursive root isolation on an interval (a, b].
 */
void sturm_seqence_isolate_roots(
    const upolynomial_dense_t* S, size_t S_size,
    algebraic_number_t* roots, size_t* roots_size,
    const dyadic_interval_t* interval,
    int a_sgn_changes, int b_sgn_changes)
{
  dyadic_interval_t I;
  dyadic_interval_ops.construct_copy(&I, interval);

  for (;;) {

    if (debug_trace_ops.is_enabled("roots")) {
      tracef("sturm_seqence_isolate_roots(");
      upolynomial_dense_ops.print(S, trace_out);
      tracef(", ");
      dyadic_interval_ops.print(&I, trace_out);
      tracef(")\n");
      tracef("a_sgn_changes = %d\n", a_sgn_changes);
      tracef("b_sgn_changes = %d\n", b_sgn_changes);
    }

    // Total roots is the difference in the sign changes
    int total_roots = a_sgn_changes - b_sgn_changes;
    assert(total_roots != 0);

    // For one root, check if we hit the root exactly
    if (total_roots == 1) {
      // Isolated one check if it's in b
      if (upolynomial_dense_ops.sgn_at_dyadic_rational(&S[0], &I.b) == 0) {
        // Copy out the interval [a, b]
        algebraic_number_ops.construct_from_dyadic_rational(&roots[*roots_size], &I.b);
        dyadic_interval_ops.destruct(&I);
        (*roots_size) ++;
        return;
      } else if (upolynomial_dense_ops.sgn_at_dyadic_rational(&S[0], &I.a) != 0) {
        // Copy out the open interval (a, b)
        I.b_open = 1;
        upolynomial_t* f = upolynomial_dense_ops.to_upolynomial(&S[0], Z);
        algebraic_number_ops.construct(&roots[*roots_size], f, &I);
        dyadic_interval_ops.destruct(&I);
        (*roots_size) ++;
        return;
      }
    }

    // Continue with the splits (a, m], (m, b]
    dyadic_interval_t I_left;
    dyadic_interval_t I_right;
    dyadic_interval_ops.construct_from_split(&I_left, &I_right, &I, 0, 1);

    // Number of sign changes in the mid-point
    int m_sgn_changes = sturm_seqence_count_sign_changes_dyadic(S, S_size, &I_left.b, a_sgn_changes);

    // No need for recursion if one of the intervals keeps all the roots
    if (a_sgn_changes == m_sgn_changes) {
      // All the roots are in the right interval
      dyadic_interval_ops.swap(&I, &I_right);
    } else if (b_sgn_changes == m_sgn_changes) {
      // All the roots are in the left interval
      dyadic_interval_ops.swap(&I, &I_left);
    } else {
      // Recurse
      // Isolate in the right interval
      sturm_seqence_isolate_roots(S, S_size, roots, roots_size, &I_left, a_sgn_changes, m_sgn_changes);
      // Isolate in the left interval
      sturm_seqence_isolate_roots(S, S_size, roots, roots_size, &I_right, m_sgn_changes, b_sgn_changes);
      // Remove the temporaries
      dyadic_interval_ops.destruct(&I);
      dyadic_interval_ops.destruct(&I_left);
      dyadic_interval_ops.destruct(&I_right);
      // Done
      return;
    }

    // Remove the temporaries
    dyadic_interval_ops.destruct(&I_left);
    dyadic_interval_ops.destruct(&I_right);
  }
}

/**
 * Same as the root counting except the we
 */
void upolynomial_roots_isolate_sturm(const upolynomial_t* f, algebraic_number_t* roots, size_t* roots_size) {

  assert(f->K == Z);

  TRACE("roots", "upolynomial_root_isolate_sturm(%P)\n", f);

  *roots_size = 0;

  // Special case for the constants
  if (upolynomial_ops.degree(f) == 0) {
    assert(!upolynomial_ops.is_zero(f));
    return;
  }

  // Get the square-free factorization and then count roots for each factor.
  upolynomial_factors_t* square_free_factors = upolynomial_ops.factor_square_free(f);

  int factor_i;
  for (factor_i = 0; factor_i < square_free_factors->size; ++ factor_i) {
    // The factor we are working with
    const upolynomial_t* factor = square_free_factors->factors[factor_i];
    int factor_deg = upolynomial_ops.degree(factor);

    // Compute the Sturm sequence for the factor
    size_t sturm_sequence_size;
    upolynomial_dense_t* sturm_sequence = malloc((factor_deg+1)*sizeof(upolynomial_dense_t));
    upolynomial_compute_sturm_sequence(factor, sturm_sequence, &sturm_sequence_size);

    // Get the total number of roots
    int total_count = sturm_seqence_count_roots_dyadic(sturm_sequence, sturm_sequence_size, 0);
    // Now grow the interval (-1, 1] until it captures all the roots
    dyadic_interval_t interval_all;
    dyadic_interval_ops.construct_from_int(&interval_all, -1, 1, 1, 1);
    int a_sgn_changes, b_sgn_changes;
    for (;;) {

      if (debug_trace_ops.is_enabled("roots")) {
        tracef("interval_all: ");
        dyadic_interval_ops.print(&interval_all, trace_out);
        tracef("\n");
      }

      // Compute the sign changes
      a_sgn_changes = sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval_all.a, sturm_sequence_size);
      b_sgn_changes = sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval_all.b, sturm_sequence_size);
      if (a_sgn_changes - b_sgn_changes == total_count) {
        break;
      }
      // Get the total number of roots
      dyadic_interval_ops.scale(&interval_all, 1);
    }

    // Isolate the roots
    if (a_sgn_changes - b_sgn_changes > 0) {
      sturm_seqence_isolate_roots(sturm_sequence, sturm_sequence_size, roots, roots_size, &interval_all, a_sgn_changes, b_sgn_changes);
    }

    // Destroy the temporaries
    dyadic_interval_ops.destruct(&interval_all);
    int i;
    for (i = 0; i < sturm_sequence_size; ++i) {
      upolynomial_dense_ops.destruct(sturm_sequence + i);
    }
    free(sturm_sequence);
  }

  TRACE("roots", "upolynomial_root_isolate_sturm(%P) = %d \n", f, *roots_size);

  // Sort the roots
  qsort(roots, *roots_size, sizeof(algebraic_number_t), algebraic_number_ops.cmp_void);

  // Destroy the factors
  upolynomial_factors_ops.destruct(square_free_factors, 1);
}
