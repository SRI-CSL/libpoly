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

#include <rational_interval.h>
#include <dyadic_interval.h>

#include "upolynomial/upolynomial.h"
#include "upolynomial/factors.h"
#include "upolynomial/root_finding.h"
#include "upolynomial/factorization.h"
#include "upolynomial/upolynomial_vector.h"
#include "upolynomial/upolynomial_dense.h"
#include "upolynomial/output.h"

#include <assert.h>
#include <stdlib.h>

#include "utils/debug_trace.h"

/** Field order up to we're doing brute force. */
#define FIELD_ORDER_LIMIT 1000

/** Negative infinity */
#define INF_N (void*) 0
/** Positive infinity */
#define INF_P (void*) 1

void upolynomial_compute_sturm_sequence(const lp_upolynomial_t* f, upolynomial_dense_t* S, size_t* size) {

  if (trace_is_enabled("roots")) {
    tracef("upolynomial_compute_sturm_sequence("); lp_upolynomial_print(f, trace_out); tracef("\n");
  }

  lp_integer_t a;
  integer_construct_from_int(lp_Z, &a, 0);

  // Min size for the polynomials
  size_t f_deg = lp_upolynomial_degree(f);

  // f[0] = pp(f)
  upolynomial_dense_construct_p(&S[0], f_deg + 1, f);
  upolynomial_dense_mk_primitive_Z(&S[0], 1);

  if (trace_is_enabled("roots")) {
    tracef("S[0] = ");
    upolynomial_dense_print(&S[0], trace_out);
    tracef("\n");
  }

  // f[1] = f'
  upolynomial_dense_construct(&S[1], f_deg + 1);
  upolynomial_dense_derivative(lp_Z, &S[0], &S[1]);
  upolynomial_dense_mk_primitive_Z(&S[1], 1);

  if (trace_is_enabled("roots")) {
    tracef("S[1] = ");
    upolynomial_dense_print(&S[1], trace_out);
    tracef("\n");
  }

  // Until we hit a constant polynomial
  int i = 1;
  while (S[i].size > 1) {
    i ++;
    upolynomial_dense_construct(&S[i], f_deg + 1);
    // Compute a*S[i-2] = div*S[i-1] + b*S[i]
    upolynomial_dense_reduce_Z(&S[i-2], &S[i-1], &a, &S[i]);
    upolynomial_dense_mk_primitive_Z(&S[i], 0);

    // If the coefficient of the reduction is not negative, we have to flip the
    // sign to get a Sturm sequence
    if (integer_sgn(lp_Z, &a) > 0) {
      upolynomial_dense_negate(&S[i], lp_Z);
    }

    if (trace_is_enabled("roots")) {
      tracef("S[%d] = ", i); upolynomial_dense_print(&S[i], trace_out); tracef("\n");
    }
  }

  // Remove the temp
  integer_destruct(&a);

  // Size
  *size = i + 1;
}

/**
 * Compute the number sgn_changes(a) given Sturm sequence and a. If a is 0 or
 * 1 as a pointer, we evaluate at -inf, +inf respectively.
 */
static
int sturm_seqence_count_sign_changes(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const lp_rational_t* a, int max_changes)
{
  int i, sgn_a, sgn_a_previous = 0, sgn_a_changes_count = 0;
  for (i = 0; i < sturm_sequence_size && sgn_a_changes_count < max_changes; ++ i) {
    // Get the sign of S[i] at a
    sgn_a = (a == INF_N) ? upolynomial_dense_sgn_at_minus_inf(&sturm_sequence[i])
          : (a == INF_P) ? upolynomial_dense_sgn_at_plus_inf(&sturm_sequence[i])
          : upolynomial_dense_sgn_at_rational(&sturm_sequence[i], a);
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
static
int sturm_seqence_count_sign_changes_dyadic(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const lp_dyadic_rational_t* a, int max_changes)
{
  int i, sgn_a, sgn_a_previous = 0, sgn_a_changes_count = 0;
  for (i = 0; i < sturm_sequence_size && sgn_a_changes_count < max_changes; ++ i) {
    // Get the sign of S[i] at a
    sgn_a = (a == INF_N) ? upolynomial_dense_sgn_at_minus_inf(&sturm_sequence[i])
          : (a == INF_P) ? upolynomial_dense_sgn_at_plus_inf(&sturm_sequence[i])
          : upolynomial_dense_sgn_at_dyadic_rational(&sturm_sequence[i], a);
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
static
int sturm_seqence_count_roots(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const lp_rational_interval_t* interval)
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
    if (interval->b_open && upolynomial_dense_sgn_at_rational(&sturm_sequence[0], &interval->b) == 0) {
      root_count --;
    }
    if (!interval->a_open && upolynomial_dense_sgn_at_rational(&sturm_sequence[0], &interval->a) != 0) {
      root_count ++;
    }
  }
  // Done
  return root_count;
}

/** Count roots in the given interval */
static
int sturm_seqence_count_roots_dyadic(
    const upolynomial_dense_t* sturm_sequence, int sturm_sequence_size,
    const lp_dyadic_interval_t* interval)
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
    if (interval->b_open && upolynomial_dense_sgn_at_dyadic_rational(&sturm_sequence[0], &interval->b) == 0) {
      root_count --;
    }
    if (!interval->a_open && upolynomial_dense_sgn_at_dyadic_rational(&sturm_sequence[0], &interval->a) != 0) {
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
int upolynomial_roots_count_sturm(const lp_upolynomial_t* f, const lp_rational_interval_t* interval) {

  assert(f->K == lp_Z);

  if (trace_is_enabled("roots")) {
    tracef("upolynomial_root_count_sturm(");
    lp_upolynomial_print(f, trace_out); tracef(", ");
    lp_rational_interval_print(interval, trace_out); tracef("\n");
  }

  int total_count = 0;

  // Special case for the constants
  if (lp_upolynomial_degree(f) == 0) {
    assert(!lp_upolynomial_is_zero(f));
    return 0;
  }

  // Get the square-free factorization and then count roots for each factor.
  lp_upolynomial_factors_t* square_free_factors = lp_upolynomial_factor_square_free(f);

  size_t factor_i;
  for (factor_i = 0; factor_i < square_free_factors->size; ++ factor_i) {
    // The factor we are working with
    const lp_upolynomial_t* factor = square_free_factors->factors[factor_i];

    // Compute the Sturm sequence for the factor
    size_t sturm_sequence_size;
    upolynomial_dense_t* sturm_sequence = malloc((lp_upolynomial_degree(factor) + 1)*sizeof(upolynomial_dense_t));
    upolynomial_compute_sturm_sequence(factor, sturm_sequence, &sturm_sequence_size);

    // Add the number of roots
    total_count += sturm_seqence_count_roots(sturm_sequence, sturm_sequence_size, interval);

    // Destroy the temporaries
    size_t i;
    for (i = 0; i < sturm_sequence_size; ++i) {
      upolynomial_dense_destruct(sturm_sequence + i);
    }
    free(sturm_sequence);
  }

  // Destroy the factors
  lp_upolynomial_factors_destruct(square_free_factors, 1);

  // Return the total number of roots
  return total_count;
}

/**
 * Recursive root isolation on an interval (a, b].
 */
static
void sturm_seqence_isolate_roots(
    const upolynomial_dense_t* S, size_t S_size,
    lp_algebraic_number_t* roots, size_t* roots_size,
    const lp_dyadic_interval_t* interval,
    int a_sgn_changes, int b_sgn_changes)
{
  lp_dyadic_interval_t I;
  lp_dyadic_interval_construct_copy(&I, interval);

  for (;;) {

    if (trace_is_enabled("roots")) {
      tracef("sturm_seqence_isolate_roots(");
      upolynomial_dense_print(S, trace_out);
      tracef(", ");
      lp_dyadic_interval_print(&I, trace_out);
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
      if (upolynomial_dense_sgn_at_dyadic_rational(&S[0], &I.b) == 0) {
        // Copy out the interval [a, b]
        lp_algebraic_number_construct_from_dyadic_rational(&roots[*roots_size], &I.b);
        lp_dyadic_interval_destruct(&I);
        (*roots_size) ++;
        return;
      } else if (upolynomial_dense_sgn_at_dyadic_rational(&S[0], &I.a) != 0) {
        // Copy out the open interval (a, b)
        I.b_open = 1;
        lp_upolynomial_t* f = upolynomial_dense_to_upolynomial(&S[0], lp_Z);
        if (trace_is_enabled("roots")) {
          tracef("f = "); lp_upolynomial_print(f, trace_out); tracef("\n");
        }
        lp_algebraic_number_construct(&roots[*roots_size], f, &I);
        lp_dyadic_interval_destruct(&I);
        (*roots_size) ++;
        return;
      }
    }

    // Continue with the splits (a, m], (m, b]
    lp_dyadic_interval_t I_left;
    lp_dyadic_interval_t I_right;
    lp_dyadic_interval_construct_from_split(&I_left, &I_right, &I, 0, 1);

    // Number of sign changes in the mid-point
    int m_sgn_changes = sturm_seqence_count_sign_changes_dyadic(S, S_size, &I_left.b, a_sgn_changes);

    // No need for recursion if one of the intervals keeps all the roots
    if (a_sgn_changes == m_sgn_changes) {
      // All the roots are in the right interval
      lp_dyadic_interval_swap(&I, &I_right);
    } else if (b_sgn_changes == m_sgn_changes) {
      // All the roots are in the left interval
      lp_dyadic_interval_swap(&I, &I_left);
    } else {
      // Recurse
      // Isolate in the right interval
      sturm_seqence_isolate_roots(S, S_size, roots, roots_size, &I_left, a_sgn_changes, m_sgn_changes);
      // Isolate in the left interval
      sturm_seqence_isolate_roots(S, S_size, roots, roots_size, &I_right, m_sgn_changes, b_sgn_changes);
      // Remove the temporaries
      lp_dyadic_interval_destruct(&I);
      lp_dyadic_interval_destruct(&I_left);
      lp_dyadic_interval_destruct(&I_right);
      // Done
      return;
    }

    // Remove the temporaries
    lp_dyadic_interval_destruct(&I_left);
    lp_dyadic_interval_destruct(&I_right);
  }
}

/**
 * Same as the root counting except the we
 */
void upolynomial_roots_isolate_sturm(const lp_upolynomial_t* f, lp_algebraic_number_t* roots, size_t* roots_size) {

  assert(f->K == lp_Z);

  if (trace_is_enabled("roots")) {
    tracef("upolynomial_root_isolate_sturm("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }

  *roots_size = 0;

  // Special case for the constants
  if (lp_upolynomial_degree(f) == 0) {
    assert(!lp_upolynomial_is_zero(f));
    return;
  }

  // Get the square-free factorization and then count roots for each factor.
  lp_upolynomial_factors_t* square_free_factors = lp_upolynomial_factor_square_free(f);

  //
  // SQUARE-FREE FACTORS CAN NOT SHARE ROOTS: we therefore us the input array directly
  //

  size_t factor_i;
  for (factor_i = 0; factor_i < square_free_factors->size; ++ factor_i) {

    // The factor we are working with
    const lp_upolynomial_t* factor = square_free_factors->factors[factor_i];
    int factor_deg = lp_upolynomial_degree(factor);

    if (trace_is_enabled("roots")) {
      tracef("upolynomial_root_isolate_sturm(): factor = "); lp_upolynomial_print(factor, trace_out); tracef(")\n");
    }

    // Check if it's a power of x
    if (!lp_upolynomial_const_term(factor)) {
      assert(factor_deg == 1);
      // Add 0 as a root
      lp_algebraic_number_construct_zero(roots + *roots_size);
      (*roots_size) ++;
      assert(*roots_size <= lp_upolynomial_degree(f));
      continue;
    }

    // Compute the Sturm sequence for the factor
    size_t sturm_sequence_size;
    upolynomial_dense_t* sturm_sequence = malloc((factor_deg+1)*sizeof(upolynomial_dense_t));
    upolynomial_compute_sturm_sequence(factor, sturm_sequence, &sturm_sequence_size);

    // Get the total number of roots
    int total_count = sturm_seqence_count_roots_dyadic(sturm_sequence, sturm_sequence_size, 0);
    // Now grow the interval (-1, 1] until it captures all the roots
    lp_dyadic_interval_t interval_all;
    lp_dyadic_interval_construct_from_int(&interval_all, -1, 1, 1, 1);
    int a_sgn_changes, b_sgn_changes;
    for (;;) {

      if (trace_is_enabled("roots")) {
        tracef("interval_all: ");
        lp_dyadic_interval_print(&interval_all, trace_out);
        tracef("\n");
      }

      // Compute the sign changes
      a_sgn_changes = sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval_all.a, sturm_sequence_size);
      b_sgn_changes = sturm_seqence_count_sign_changes_dyadic(sturm_sequence, sturm_sequence_size, &interval_all.b, sturm_sequence_size);
      if (a_sgn_changes - b_sgn_changes == total_count) {
        break;
      }
      // Get the total number of roots
      lp_dyadic_interval_scale(&interval_all, 1);
    }

    // Isolate the roots
    if (a_sgn_changes - b_sgn_changes > 0) {
      size_t current_roots_size = 0;
      sturm_seqence_isolate_roots(sturm_sequence, sturm_sequence_size, roots + *roots_size, &current_roots_size, &interval_all, a_sgn_changes, b_sgn_changes);
      (*roots_size) += current_roots_size;
      assert(*roots_size <= lp_upolynomial_degree(f));
    }

    // Destroy the temporaries
    lp_dyadic_interval_destruct(&interval_all);
    for (size_t i = 0; i < sturm_sequence_size; ++i) {
      upolynomial_dense_destruct(sturm_sequence + i);
    }
    free(sturm_sequence);
  }

  if (trace_is_enabled("roots")) {
    tracef("upolynomial_root_isolate_sturm(");
    lp_upolynomial_print(f, trace_out);
    tracef(" = %zu \n", *roots_size);
  }

  // Sort the roots
  qsort(roots, *roots_size, sizeof(lp_algebraic_number_t), lp_algebraic_number_cmp_void);

  // Destroy the factors
  lp_upolynomial_factors_destruct(square_free_factors, 1);
}

/**
 * helper function for rabin root finding
 * returns x^p mod f
 */
static
lp_upolynomial_t* upolynomial_power_mod(const lp_upolynomial_t *b, const lp_integer_t *e, const lp_upolynomial_t *m) {
  lp_integer_t ee;
  lp_integer_construct_copy(lp_Z, &ee, e);
  lp_upolynomial_t *acc = lp_upolynomial_construct_power(m->K, 0, 1);
  lp_upolynomial_t *power = lp_upolynomial_construct_copy(b);

  while (!lp_integer_is_zero(lp_Z, &ee)) {
    if (integer_is_odd(&ee)) {
      upolynomial_op_inplace(lp_upolynomial_mul, &acc, power);
      upolynomial_op_inplace(lp_upolynomial_rem_exact, &acc, m);
    }
    upolynomial_op_inplace(lp_upolynomial_mul, &power, power);
    upolynomial_op_inplace(lp_upolynomial_rem_exact, &power, m);

    integer_div_floor_pow2(&ee, &ee, 1);
  }

  lp_integer_destruct(&ee);
  lp_upolynomial_delete(power);

  return acc;
}

/**
 * helper function for rabin root finding
 * returns gcd(f, x^p - x)
 */
static
lp_upolynomial_t* upolynomial_reduce_to_linear(const lp_upolynomial_t *f ) {
  lp_upolynomial_t *x = lp_upolynomial_construct_power(f->K, 1, 1);
  lp_upolynomial_t *pm = upolynomial_power_mod(x, &f->K->M, f);
  lp_upolynomial_t *fp = lp_upolynomial_sub(pm, x);
  lp_upolynomial_t *linear = lp_upolynomial_gcd(f, fp);

  lp_upolynomial_delete(x);
  lp_upolynomial_delete(pm);
  lp_upolynomial_delete(fp);

  return linear;
}

/**
 * Implements Rabin root-finding.
 * Rabin, Michael O. "Probabilistic algorithms in finite fields." SIAM Journal on computing 9.2 (1980): 273-280.
 */
static
void upolynomial_roots_find_rabin(const lp_upolynomial_t* f, lp_integer_t* roots, size_t* roots_size) {
  const lp_int_ring_t* K = f->K;

  *roots_size = 0;

  gmp_randstate_t state;
  gmp_randinit_default(state);

  lp_integer_t s;
  lp_integer_construct_copy(lp_Z, &s, &K->M);
  // right shift by p1
  integer_div_floor_pow2(&s, &s, 1);

  // construct helper polynomials (1 and x)
  lp_upolynomial_t *p1 = lp_upolynomial_construct_power(K, 0, 1);
  lp_upolynomial_t *px = lp_upolynomial_construct_power(K, 1, 1);

  lp_upolynomial_vector_t *to_factor = lp_upolynomial_vector_construct();
  lp_upolynomial_vector_move_back(to_factor, upolynomial_reduce_to_linear(f));

  while (lp_upolynomial_vector_size(to_factor)) {
    lp_upolynomial_t *p = lp_upolynomial_vector_pop(to_factor);
    size_t d = lp_upolynomial_degree(p);
    bool const_is_zero = !lp_upolynomial_const_term(p) || lp_integer_is_zero(lp_Z, lp_upolynomial_const_term(p));
    if (d == 0) {
      // it's dead
    } else if (const_is_zero) {
      // it has a zero root
      assert(lp_upolynomial_degree(f) > *roots_size);
      lp_integer_construct_from_int(lp_Z, &roots[(*roots_size)++], 0);
      lp_upolynomial_vector_move_back(to_factor, lp_upolynomial_div_exact(p, px));
    } else if (d == 1) {
      // push the zero
      lp_upolynomial_make_monic_in_place(p);
      lp_upolynomial_neg_in_place(p);
      assert(lp_upolynomial_const_term(p));
      assert(lp_upolynomial_degree(f) > *roots_size);
      lp_integer_construct_copy(lp_Z, &roots[(*roots_size)++], lp_upolynomial_const_term(p));
    } else {
      // gen (px - delta)
      lp_upolynomial_t *xd = lp_upolynomial_construct_from_int(K, 1, (int[2]){1, 1});
      lp_integer_t *delta = (lp_integer_t*) lp_upolynomial_const_term(xd);
      assert(delta);

      bool done = false;
      while (!done) {
        // set constant term delta in xd to random value
        mpz_urandomm(delta, state, &K->M);
        integer_ring_normalize(K, delta);
        integer_neg(K, delta, delta);
        // calculate gcd
        lp_upolynomial_t *h = upolynomial_power_mod(xd, &s, p);
        upolynomial_op_inplace(lp_upolynomial_sub, &h, p1);
        upolynomial_op_inplace(lp_upolynomial_gcd, &h, p);
        if (0 < lp_upolynomial_degree(h) && lp_upolynomial_degree(h) < d) {
          lp_upolynomial_vector_push_back(to_factor, h);
          lp_upolynomial_vector_move_back(to_factor, lp_upolynomial_div_exact(p, h));
          done = true;
        }
        lp_upolynomial_delete(h);
      }
      lp_upolynomial_delete(xd);
    }
    lp_upolynomial_delete(p);
  }

  lp_integer_destruct(&s);
  lp_upolynomial_delete(p1);
  lp_upolynomial_delete(px);
  lp_upolynomial_vector_delete(to_factor);
  gmp_randclear(state);
}

/**
 * iterates over all p possible values of Zp and check each
 */
static
void upolynomial_roots_find_brute_force(const lp_upolynomial_t* f, lp_integer_t* roots, size_t* roots_size) {
  const lp_int_ring_t* K = f->K;

  size_t d = lp_upolynomial_degree(f);
  assert(d > 0);

  assert(lp_integer_fits_int(&K->M));
  long p = integer_to_int(&K->M);
  assert(p < FIELD_ORDER_LIMIT);

  *roots_size = 0;

  // Values
  lp_integer_t x, value;
  integer_construct_copy(lp_Z, &x, &K->lb);
  integer_construct(&value);

  for (long x_int = 0; x_int < p; ++x_int) {
    assert(integer_in_ring(K, &x));

    // Evaluate the polynomial
    lp_upolynomial_evaluate_at_integer(f, &x, &value);

    // If zero we found one
    if (integer_sgn(K, &value) == 0) {
      lp_integer_construct_copy(K, roots + *roots_size, &x);
      (*roots_size) ++;

      // check if we found all
      if (*roots_size == d) {
        break;
      }
    }

    // next element
    integer_inc(lp_Z, &x);
  }
  assert(*roots_size == d || lp_integer_cmp(lp_Z, &x, &K->ub));

  integer_destruct(&value);
  integer_destruct(&x);
}

void upolynomial_roots_find_Zp(const lp_upolynomial_t* f, lp_integer_t** roots, size_t* roots_size) {
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_roots_find_Zp("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }

  const lp_int_ring_t* K = f->K;
  assert(K && K->is_prime);

  size_t d = lp_upolynomial_degree(f);

  *roots = malloc(sizeof(lp_integer_t) * d);

  // depending on the finite field size, choose appropriate function
  if (integer_cmp_int(lp_Z, &K->M, FIELD_ORDER_LIMIT) < 0) {
    upolynomial_roots_find_brute_force(f, *roots, roots_size);
  } else {
    upolynomial_roots_find_rabin(f, *roots, roots_size);
  }

  // realloc has undefined behavior in case size is zero
  if (*roots_size == 0) {
    free(*roots);
    *roots = NULL;
  } else if (*roots_size < d) {
    *roots = realloc(*roots, sizeof(lp_integer_t) * (*roots_size));
  }

  if (trace_is_enabled("roots")) {
    tracef("upolynomial_roots_find_Zp() found %zu roots\n", *roots_size);
  }
}
