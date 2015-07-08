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

#include "upolynomial/factorization.h"
#include "upolynomial/upolynomial.h"
#include "upolynomial/factors.h"
#include "upolynomial/bounds.h"
#include "upolynomial/output.h"

#include "utils/statistics.h"
#include "utils/debug_trace.h"

#include <assert.h>
#include <stdlib.h>


STAT_DECLARE(int, upolynomial, factor_square_free)
STAT_DECLARE(int, upolynomial, factor_distinct_degree)
STAT_DECLARE(int, upolynomial, factor_berlekamp_square_free)

/**
 * We are given a polynomial f and we will return it's square-free factorization
 *
 *  f = \prod_{k} f_k^k
 *
 * where each f_k is square free and all gcd(f_i, f_j) = 1.
 *
 * First, note that if f has a square factor, then the gcd(f, f') != 1. This is
 * because if f = u^2*v, then f' = 2*u*u'*v + u^2*v, so gcd(f, f') would include
 * at least u != 1.
 *
 * Now,consider f'
 *
 *  f' = \sum_{k} k*f_k^{k-1}*f_k'*\prod_{i!=k} f_i^i (if we are in Z_p)
 *     = \sum_{p \ndiv k} k*f_k^{k-1}*f_k'*\prod_{i!=k}
 *
 * It is easy to see that f' is divisible by each f_k^{k-1} for k such that
 * p \ndiv k, but not f_k^k. It is also divisible by f_k^k for k such that p
 * \div k. Since these are the only possible factors of f, then
 *
 *  P = gcd(f, f') = \prod_{p \ndiv k} f_k^{k-1} * \prod_{p \div k} f_k^k
 *
 * With P as above, then
 *
 *  L = f/P = \prod_{p \ndiv k} f_k               -> Linear product (can be 1)
 *
 * To compute the first factor we then compute (loop start)
 *
 *  R = gcd(P, L) = \prod_{k > 1 & p \ndiv k} f_k -> Rest of linear product
 *  O = L / R                                     -> First factor (output)
 *
 * To continue the loop, we set
 *
 *  P = P / R     -> Reduced powers of f_k with p \ndiv k, k > 1
 *  L = R         -> Linear product of f_k with p \ndiv k, k > 1
 *
 * And go on with the loop.
 *
 * The loop ends when L = 1 and then we know that P = \prod_{p \div k} f_k^k
 * which we know how to special case (if P != 1).
 */
lp_upolynomial_factors_t* lp_upolynomial_factor_square_free_primitive(const lp_upolynomial_t* f) {

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_square_free("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }
  STAT(upolynomial, factor_square_free) ++;

  assert(!f->K || !f->K->is_prime || lp_upolynomial_is_monic(f));
  assert(f->K || lp_upolynomial_is_primitive(f));
  assert(lp_upolynomial_const_term(f));

  lp_upolynomial_factors_t* factors = 0;

  // Special case for constants
  if (lp_upolynomial_degree(f) == 0) {
    const lp_integer_t* c = lp_upolynomial_const_term(f);
    assert(c);
    factors = lp_upolynomial_factors_construct();
    integer_assign(f->K, &factors->constant, c);
    return factors;
  }

  // Derivative
  lp_upolynomial_t* d_f = lp_upolynomial_derivative(f);

  if (lp_upolynomial_is_zero(d_f)) {
    assert(f->K && f->K->is_prime);
    // f' is zero for a non-zero polynomial => f has to be of the form
    // f = \sum a_k x^(p*d_k) = f_p(x^p) where f_p = \sum a_k x^d_k
    // we factor f_p and then return f_p(x^p)=(f_p)^p
    int p = integer_to_int(&f->K->M);
    lp_upolynomial_t* f_p = lp_upolynomial_div_degrees(f, p);
    factors = lp_upolynomial_factor_square_free_primitive(f_p);
    size_t i;
    for (i = 0; i < factors->size; ++ i) {
      factors->multiplicities[i] *= p;
    }
    lp_upolynomial_delete(f_p);
  } else {

    // Construct the factorization
    factors = lp_upolynomial_factors_construct();

    // Degree of the factor
    int k = 1;
    // P = GCD(f, f')
    lp_upolynomial_t* P = lp_upolynomial_gcd(f, d_f);
    if (trace_is_enabled("factorization")) {
      tracef("P = "); lp_upolynomial_print(P, trace_out); tracef("\n");
    }
    // L = f/P
    lp_upolynomial_t* L = lp_upolynomial_div_exact(f, P);
    if (trace_is_enabled("factorization")) {
      tracef("L = "); lp_upolynomial_print(L, trace_out); tracef("\n");
    }

    while (lp_upolynomial_degree(L) > 0) {
      // R = gcd(P, L)
      lp_upolynomial_t* R = lp_upolynomial_gcd(P, L);
      if (trace_is_enabled("factorization")) {
        tracef("R = "); lp_upolynomial_print(R, trace_out); tracef("\n");
      }
      // O = L / R (it can be constant if there is no factor of power k)
      if (lp_upolynomial_cmp(L, R)) {
        lp_upolynomial_t* O = lp_upolynomial_div_exact(L, R);
        if (trace_is_enabled("factorization")) {
          tracef("O = "); lp_upolynomial_print(O, trace_out); tracef("\n");
        }
        // Record the output
        lp_upolynomial_factors_add(factors, O, k);
      }
      // P = P / R
      lp_upolynomial_t* tmp = P;
      P = lp_upolynomial_div_exact(P, R);
      if (trace_is_enabled("factorization")) {
        tracef("P = "); lp_upolynomial_print(P, trace_out); tracef("\n");
      }
      lp_upolynomial_delete(tmp);
      // L = R
      lp_upolynomial_delete(L);
      L = R;
      if (trace_is_enabled("factorization")) {
        tracef("L = "); lp_upolynomial_print(L, trace_out); tracef("\n");
      }
      // Next degree
      k = k + 1;
    }

    // If P has content, it is a power of p
    if (lp_upolynomial_degree(P) > 0) {
      int p = integer_to_int(&f->K->M);
      lp_upolynomial_t* P_p = lp_upolynomial_div_degrees(P, p);
      lp_upolynomial_factors_t* sub_factors = lp_upolynomial_factor_square_free_primitive(P_p);
      size_t i;
      for (i = 0; i < sub_factors->size; ++ i) {
        lp_upolynomial_factors_add(factors, sub_factors->factors[i], sub_factors->multiplicities[i] * p);
      }
      lp_upolynomial_factors_destruct(sub_factors, 0);
      lp_upolynomial_delete(P_p);
    }

    lp_upolynomial_delete(P);
    lp_upolynomial_delete(L);
  }

  // Destroy
  lp_upolynomial_delete(d_f);

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_square_free("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  return factors;
}

/**
 * Make sure its primitive and then discharge.
 */
lp_upolynomial_factors_t* lp_upolynomial_factor_square_free(const lp_upolynomial_t* f) {

  lp_integer_t content;
  lp_integer_construct(&content);

  // We factor the primitive part
  lp_upolynomial_t* f_pp = 0;

  if (f->K == lp_Z) {
    // The content of the polynomial (constant of the factorization)
    lp_upolynomial_content_Z(f, &content);
    f_pp = lp_upolynomial_primitive_part_Z(f);
  } else {
    assert(f->K->is_prime);
    // Make monic
    integer_assign(lp_Z, &content, lp_upolynomial_lead_coeff(f));
    f_pp = lp_upolynomial_div_exact_c(f, &content);
  }

  lp_upolynomial_factors_t* sq_free_factors = 0;

  // Take out the power of x^k
  if (lp_upolynomial_const_term(f_pp)) {
    // Get a square-free decomposition of f
    sq_free_factors = lp_upolynomial_factor_square_free_primitive(f_pp);
  } else {
    // Get a copy without the power of x
    lp_upolynomial_t* f_pp_nonzero = lp_upolynomial_construct_copy(f_pp);
    size_t x_degree = f_pp_nonzero->monomials[0].degree;
    size_t i = 0;
    for (i = 0; i < f_pp_nonzero->size; ++ i) {
      f_pp_nonzero->monomials[i].degree -= x_degree;
    }
    // Get a square-free decomposition of f
    sq_free_factors = lp_upolynomial_factor_square_free_primitive(f_pp_nonzero);
    // Add x^k to the factorization
    lp_upolynomial_t* x_upoly = lp_upolynomial_construct_power(f->K, 1, 1);
    lp_upolynomial_factors_add(sq_free_factors, x_upoly, x_degree);
    // Delete the temp
    lp_upolynomial_delete(f_pp_nonzero);
  }

  // Add the constant
  integer_mul(f->K, &sq_free_factors->constant, &sq_free_factors->constant, &content);

  // Remove the temps
  integer_destruct(&content);
  lp_upolynomial_delete(f_pp);

  // Return the result
  return sq_free_factors;
}


/**
 * We are given a monic, square-free polynomial f in Z_p and we will return
 * its distinct degree factorization
 *
 *  f = f_d_1 * f_d_2 * ... f_d_k
 *
 * where d_1 < d_2 < ... < d_k, and each f_d_i is a product of irreducible
 * polynomials of degree d_i.
 *
 * Theorem: Let I(d) be the set of all irreducible polynomials g in Z_p[x] with
 * deg(g) = d, then
 *
 *                ----
 *  x^(p^m) - x = |  | g(x)
 *                |  |
 *             g \in I(m)
 *                d|m
 *
 * We can therefore enumerate d and get the distinct d-factors by filtering
 * with gcd through  x^(p^d)-x.
 */
lp_upolynomial_factors_t* upolynomial_factor_distinct_degree(const lp_upolynomial_t* f) {

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_distinct_degree("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }
  STAT(upolynomial, factor_distinct_degree) ++;

  lp_int_ring_t* K = f->K;
  assert(K && K->is_prime);
  assert(lp_upolynomial_is_monic(f));

  // The prime
  int p = integer_to_int(&K->M);
  assert(p < 100);

  lp_upolynomial_factors_t* factors = lp_upolynomial_factors_construct();

  // Enumerate with d
  size_t d = 0;
  // Keep x TODO: optimize when switch to dense representation
  lp_upolynomial_t* x = lp_upolynomial_construct_power(K, 1, 1);
  // Keep x^p^d, start with x^p^0
  lp_upolynomial_t* x_pow = lp_upolynomial_construct_power(K, 1, 1);
  // Factors of f with deg >= d
  lp_upolynomial_t* f_rest = lp_upolynomial_construct_copy(f);

  // Temps
  lp_upolynomial_t* tmp = 0;

  do {

    // Our current degree left
    size_t f_rest_deg = lp_upolynomial_degree(f_rest);

    // If left with trivial or no two factors with deg >= d will fit, we're done
    if (f_rest_deg == 0 || 2*d > f_rest_deg) {
      if (f_rest_deg > 0) {
        lp_upolynomial_factors_add(factors, f_rest, d);
      }
      break;
    }

    // Go on to the next one
    d = d + 1;
    // Power up TODO: optimize power, make one modular
    tmp = x_pow;
    x_pow = lp_upolynomial_pow(x_pow, p);
    lp_upolynomial_delete(tmp);
    tmp = x_pow;
    x_pow = lp_upolynomial_rem_exact(x_pow, f_rest);
    lp_upolynomial_delete(tmp);

    // Compute x^q - x (big product from description)
    tmp = lp_upolynomial_sub(x_pow, x);

    // Filter with gcd
    lp_upolynomial_t* f_d = lp_upolynomial_gcd(tmp, f_rest);
    lp_upolynomial_delete(tmp);

    if (lp_upolynomial_degree(f_d) > 0) {
      // Remove the factor
      tmp = f_rest;
      f_rest = lp_upolynomial_div_exact(f_rest, tmp);
      lp_upolynomial_delete(tmp);
      // Simplify the power
      tmp = x_pow;
      x_pow = lp_upolynomial_rem_exact(x_pow, f_rest);
      lp_upolynomial_delete(tmp);
      // Remember the factor
      lp_upolynomial_factors_add(factors, f_d, d);
    }

    lp_upolynomial_delete(f_d);

  } while (1);

  // Remove temps
  lp_upolynomial_delete(f_rest);
  lp_upolynomial_delete(x);
  lp_upolynomial_delete(x_pow);
  lp_upolynomial_delete(tmp);

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_distinct_degree("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  return factors;
}


static void Q_construct(lp_integer_t* Q, size_t size, const lp_upolynomial_t* u) {

  lp_int_ring_t* K = lp_upolynomial_ring(u);
  size_t p = (size_t) integer_to_int(&K->M);

  size_t k;

  for (k = 0; k < size*size; ++ k) {
    integer_construct_from_int(lp_Z, Q + k, 0);
  }

  // Q
  for (k = 0; k < size; ++ k) {
    // x^{pk} congruent with q_{k, n-1}x^{n-1} + ... + q_{k, 0} (mod u)
    lp_upolynomial_t* x_pk = lp_upolynomial_construct_power(K, k*p, 1);
    lp_upolynomial_t* q_k = lp_upolynomial_rem_exact(x_pk, u);
    lp_upolynomial_unpack(q_k, Q + k*size);
    lp_upolynomial_delete(q_k);
    lp_upolynomial_delete(x_pk);
  }

  // Q - I
  lp_integer_t tmp, one;
  integer_construct_from_int(K, &tmp, 0);
  integer_construct_from_int(K, &one, 1);
  for(k = 0; k < size; ++ k) {
    integer_sub(K, &tmp, Q + k*size + k, &one);
    integer_swap(&tmp, Q + k*size + k);
  }
  integer_destruct(&tmp); integer_destruct(&one);
}

static void Q_column_multiply(lp_int_ring_t* K, lp_integer_t* Q, size_t size, int j, const lp_integer_t* m) {

  if (trace_is_enabled("nullspace")) {
    tracef("Multiplying column %d of Q with ", j); integer_print(m, trace_out); tracef("\n");
    tracef("Q = \n");
    integer_print_matrix(Q, size, size, trace_out);
  }

  lp_integer_t tmp;
  integer_construct_from_int(K, &tmp, 0);

  size_t k;
  for (k = 0; k < size; ++ k) {
    integer_mul(K, &tmp, Q + k*size + j, m);
    integer_swap(&tmp, Q + k*size + j);
  }

  integer_destruct(&tmp);

  TRACE("nullspace", "Q = \n");
  TRACE_CMD("nullspace", integer_print_matrix(Q, size, size, trace_out));
}

// add column j into i multiplied by m
static void Q_column_add(lp_int_ring_t* K, lp_integer_t* Q, size_t size, int i, const lp_integer_t* m, int j) {

  assert(i != j);

  // Nothing to multiply with
  if (integer_sgn(K, m) == 0) {
    return;
  }

  lp_integer_t m_copy, mul, add;

  integer_construct_copy(lp_Z, &m_copy, m);
  integer_construct_from_int(K, &mul, 0);
  integer_construct_from_int(K, &add, 0);

  if (trace_is_enabled("nullspace")) {
    tracef("Adding [%d]*", j); integer_print(m, trace_out); tracef(" of Q to %i\n", i);
    tracef("Q = \n"); integer_print_matrix(Q, size, size, trace_out);
  }

  size_t k;
  for (k = 0; k < size; ++ k) {
    integer_mul(K, &mul, Q + k*size + j, &m_copy);
    integer_add(K, &add, Q + k*size + i, &mul);
    integer_swap(&add, Q + k*size + i);
  }

  integer_destruct(&add);
  integer_destruct(&mul);
  integer_destruct(&m_copy);

  if (trace_is_enabled("nullspace")) {
    tracef("Q = \n"); integer_print_matrix(Q, size, size, trace_out);
  }
}

static void Q_null_space(lp_int_ring_t* K, lp_integer_t* Q, size_t size, lp_integer_t** v, size_t* v_size) {

  *v_size = 0;

  size_t i, j, k;

  int c[size]; // per column: which row contains the non-eliminated
  for (k = 0; k < size; ++ k) { c[k] = -1; }

  // Construct the vectors representing the null-space of Q, i.e.
  // vectors v such that v*Q = (0, ..., 0).

  // Step 1. Multiplying a column doesn't change the null-space or vectors.
  // Same holds for adding a column to another column. We therefore make it
  // "triangular"
  int pivot_found;
  lp_integer_t tmp1, tmp2;
  integer_construct_from_int(K, &tmp1, 0);
  integer_construct_from_int(K, &tmp2, 0);
  for (k = 0; k < size; ++ k) {
    pivot_found = 0;
    for (j = 0; j < size; ++ j) {
      const lp_integer_t* a_kj = Q + k*size + j;
      if (c[j] == -1 && integer_sgn(K, a_kj)) {
        // Multiply column j with the inverse to get -1
        integer_inv(K, &tmp1, a_kj);
        integer_neg(K, &tmp2, &tmp1);
        Q_column_multiply(K, Q, size, j, &tmp2);
        c[j] = k;
        // Eliminate the row
        for (i = 0; i < size; ++ i) {
          if (i != j) {
            Q_column_add(K, Q, size, i, Q + k*size + i, j);
          }
        }
        // Found the pivot
        pivot_found = 1;
      }
    }
    if (!pivot_found) {
      // Construct the vector (v0, ..., v_s) Q = (0, ..., 0):
      // * For j != k with c_j != -1, we pick Q_kj
      // * For j = k (c_k == -1) we pick 1
      // * Otherwise pick 0
      lp_integer_t* current = malloc(sizeof(lp_integer_t)*size);
      v[*v_size] = current;
      for (j = 0; j < size; ++ j) {
        if (j == k) {
          integer_construct_from_int(lp_Z, current + j, 1);
        } else {
          integer_construct_from_int(lp_Z, current + j, 0);
        }
      }
      for (j = 0; j < size; ++ j) {
        if (c[j] != -1) {
          integer_assign(lp_Z, current + c[j], Q + k*size + j);
        }
      }
      (*v_size) ++;
    }
  }

  integer_destruct(&tmp1);
  integer_destruct(&tmp2);

}

static void Q_destruct(lp_integer_t* Q, size_t size) {
  size_t i;
  for (i = 0; i < size*size; ++ i) {
    integer_destruct(Q + i);
  }
}

/**
 * Factors a given monic square-free polynomial f in Z_p using Berlekamp's
 * algorithm. F is square free and so it is of the form f = f_1*f_2*...*f_r,
 * with each f_i monic, irreducible and hence gcd(f_i, f_j) = 1.
 *
 * Let Q be the matrix:
 *
 *  | q_0,0    q_0,1    ...    q_0,n-1   |
 *  |   .        .                .      |
 *  |   .        .      ...       .      |
 *  |   .        .                .      |
 *  | q_n-1,0  q_n-1,1  ...    q_n-1,n-1 |
 *
 * where
 *
 * [1] x^pk = q_k,n-1 * x^(n-1) + ... + q_k,0 (mod f)
 *
 * Now, consider a polynomial v = v_n-1 * x^(n-1) + ... v_0 such that
 *
 * [2] (v_0, ..., v_n-1)*Q = (v_0, ..., v_n-1)
 *
 * In other words, v is in the nullspace of (Q-I). Above holds iff
 *
 * [3] v = \sum v_j x^j = \sum_j \sum_k (v_k*q_k,j) = \sum_k v_k (\sum_j q_k,j x^j)
 *          [1] (mod f) = \sum v_k x^pk = v(x^p) = v(x)^p
 *
 * Therefore for such v we have that
 *
 * [4] (v^p - v) = 0 (mod f).
 *
 * On the other hand, we know that
 *
 * [5] (v^p - v) = (v - 0)*(v - 1)* ... *(v-(p-1) (mod p)
 *
 * Combining [4,5] we get that
 *
 * [6] (v - 0)*(v - 1)* ... *(v-(p-1)) = 0 (mod f)
 *
 * Above [6] holds iff each f_i divides one of the (v - s_i) of [6]. Therefore
 * we have
 *
 *  f = \prod gcd(v - s, f), for 0 <= s < p.
 *
 * To get the factorization we find the null-space basis of Q, enumerate
 * the vectors v, and use v to split f into factors.
 *
 * Each factor gcd(v - s, f) might have further factors, which we can factor
 * further using a different null-space vector. The factors will eventually be
 * extracted since for s_i != s_j some (v - s_i) will be divisible by f_i but
 * not by f_j. In essence, each v filters the factorization further.
 *
 * Note that the first row of Q is (1, 0, ..., 0) and therefore the first
 * null-space basis vector of (Q - I) will be (1, 0, ..., 0). This means that
 * the basis size is at least 1, and the first vector (corresponding to v(x) = 1
 * does not filter the factorization, so we start at 2. If the size of the
 * basis is 1, then the polynomial is already irreducible.
 */
lp_upolynomial_factors_t* upolynomial_factor_berlekamp_square_free(const lp_upolynomial_t* f) {

  if (trace_is_enabled("berlekamp")) {
    tracef("upolynomial_factor_berlekamp_square_free("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }
  STAT(upolynomial, factor_berlekamp_square_free) ++;

  lp_upolynomial_factors_t* factors = lp_upolynomial_factors_construct();

  // Polynomial to factor and it's degree
  size_t deg_f = lp_upolynomial_degree(f);
  lp_int_ring_t* K = lp_upolynomial_ring(f);

  // If degree < 2 we're done
  if (deg_f < 2) {
    lp_upolynomial_factors_add(factors, lp_upolynomial_construct_copy(f), 1);
  } else {

    // The prime (should be small)
    assert(integer_cmp_int(lp_Z, &K->M, 100) < 0);
    int p = integer_to_int(&K->M);

    // Construct the Q matrix
    lp_integer_t Q[deg_f * deg_f];
    Q_construct(Q, deg_f, f);

    TRACE("nullspace", "Q = \n");
    TRACE_CMD("berlekamp", integer_print_matrix(Q, deg_f, deg_f, trace_out));

    // Compute the null space of Q
    lp_integer_t* v[deg_f]; // Vectors v_i
    size_t v_size = 0; // Size of the basis (and number of factors)
    Q_null_space(K, Q, deg_f, v, &v_size);

    if (trace_is_enabled("nullspace")) {
      size_t i;
      for (i = 0; i < v_size; ++i) {
        tracef("v[%zu] = \n", i);
        integer_print_matrix(v[i], 1, deg_f, trace_out);
      }
    }

    // Start with the trivial factorization
    lp_upolynomial_factors_add(factors, lp_upolynomial_construct_copy(f), 1);

    // Filter through the null-basis
    size_t i;
    int done = 0;
    for (i = 1; !done && i < v_size; ++i) {

      // v-s polynomial
      lp_upolynomial_t* v_poly[p];

      // Construct the v[i](x) - s polynomials
      int s;
      for (s = 0; s < p; ++ s) {
        v_poly[s] = lp_upolynomial_construct(K, deg_f - 1, v[i]);
        integer_dec(K, &v[i][0]);
      }

      int factor_i, factor_i_end;
      for (factor_i = 0, factor_i_end = factors->size; !done && factor_i < factor_i_end; ++ factor_i) {

        // Current factor to filter
        lp_upolynomial_t* factor = factors->factors[factor_i];

        // Maybe we're done already
        if (lp_upolynomial_degree(factor) < 2) {
          continue;
        }

        // Filter the current factorization with v_i
        int s;
        for (s = 0; !done && s < p; ++ s) {
          // Compute the gcd
          lp_upolynomial_t* gcd = lp_upolynomial_gcd(factor, v_poly[s]);
          // We only split if gcd splits factor
          size_t gcd_degree = lp_upolynomial_degree(gcd);
          if (gcd_degree > 0 && gcd_degree != lp_upolynomial_degree(factor)) {
            // Split the current factor into factor/gcd and gcd (added to back)
            lp_upolynomial_t* reduced_factor = lp_upolynomial_div_exact(factor, gcd);
            lp_upolynomial_delete(factor);
            factor = factors->factors[factor_i] = reduced_factor;
            lp_upolynomial_factors_add(factors, gcd, 1);
            // Are we done
            done = factors->size == v_size;
          } else {
            lp_upolynomial_delete(gcd);
          }
        }
      }

      // Free  the v[i](x) -s polynomials
      for (s = 0; s < p; ++ s) {
        lp_upolynomial_delete(v_poly[s]);
      }
    }

    // Deallocate the null vectors
    size_t j;
    for (i = 0; i < v_size; ++i) {
      for (j = 0; j < deg_f; ++j) {
        integer_destruct(v[i] + j);
      }
      free(v[i]);
    }

    TRACE("berlekamp", "Q = \n"); TRACE_CMD("berlekamp", integer_print_matrix(Q, deg_f, deg_f, trace_out));

    assert(v_size == factors->size);

    // Free the temp stuff
    Q_destruct(Q, deg_f);
  }

  if (trace_is_enabled("berlekamp")) {
    tracef("upolynomial_factor_berlekamp_square_free("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  // Return the result
  return factors;
}

/**
 * Factors the given polynomial f using Berlekamp's algorithm. Berlekamp's
 * algorithms works in Z_p for a prime p, and works over square-free
 * polynomials. The function first performs the square-free factorization
 * of f, and then uses Berlekamp's algorithm to factor each of the square-free
 * factors.
 */
lp_upolynomial_factors_t* upolynomial_factor_Zp(const lp_upolynomial_t* f) {

  if (trace_is_enabled("berlekamp")) {
    tracef("upolynomial_factor_Zp("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }

  lp_int_ring_t* K = f->K;

  assert(K && K->is_prime);
  assert(lp_upolynomial_degree(f) > 0);

  // Our result factorization
  lp_upolynomial_factors_t* result = lp_upolynomial_factors_construct();

  lp_upolynomial_t* to_factor = 0;
  if (lp_upolynomial_is_monic(f)) {
    to_factor = lp_upolynomial_construct_copy(f);
  } else {
    integer_assign(lp_Z, &result->constant, lp_upolynomial_lead_coeff(f));
    to_factor = lp_upolynomial_div_exact_c(f, &result->constant);
  }

  // Compute the square-free factors of f
  lp_upolynomial_factors_t* sq_free_factors = lp_upolynomial_factor_square_free(to_factor);

  // Go through the square-free factors break them apart
  int i, i_end;
  for (i = 0, i_end = sq_free_factors->size; i < i_end; ++ i) {

    // The square-free factor and its multiplicity
    lp_upolynomial_t* f_i = sq_free_factors->factors[i];
    size_t f_i_multiplicity = sq_free_factors->multiplicities[i];

    // Extract linear factors
    int x_int, p = integer_to_int(&K->M);
    lp_upolynomial_t* linear_factors_product = 0;
    for (x_int = 0; x_int < p; ++ x_int) {
      // Values
      lp_integer_t value, x;
      integer_construct_from_int(K, &value, 0);
      integer_construct_from_int(K, &x, x_int);

      // Evaluate the polynomial
      lp_upolynomial_evaluate_at_integer(f_i, &x, &value);

      // If zero we have a factor
      if (integer_sgn(lp_Z, &value) == 0) {
        int coeff[2] = { 0, 1 };
        coeff[0] = -x_int;
        lp_upolynomial_t* linear_factor = lp_upolynomial_construct_from_int(K, 1, coeff);
        lp_upolynomial_factors_add(result, linear_factor, f_i_multiplicity);
        if (linear_factors_product) {
          lp_upolynomial_t* tmp = linear_factors_product;
          linear_factors_product = lp_upolynomial_mul(linear_factors_product, linear_factor);
          lp_upolynomial_delete(tmp);
        } else {
          linear_factors_product = lp_upolynomial_construct_copy(linear_factor);
        }
      }

      integer_destruct(&x);
      integer_destruct(&value);
    }

    // Remove any linear factors
    if (linear_factors_product) {
      f_i = lp_upolynomial_div_exact(f_i, linear_factors_product);
      lp_upolynomial_delete(sq_free_factors->factors[i]);
      sq_free_factors->factors[i] = f_i;
      lp_upolynomial_delete(linear_factors_product);
    }

    if (!lp_upolynomial_is_one(f_i)) {
      // Factor it
      lp_upolynomial_factors_t* f_i_factors = upolynomial_factor_berlekamp_square_free(f_i);
      // Copy the factorization (all monic, no constant to worry about)
      size_t k;
      for (k = 0; k < f_i_factors->size; ++k) {
        assert(f_i_factors->multiplicities[k] == 1);
        lp_upolynomial_factors_add(result, f_i_factors->factors[k], f_i_multiplicity);
      }
      // Get rid of the factorization
      lp_upolynomial_factors_destruct(f_i_factors, 0);
    }
  }

  // Free the temporaries
  lp_upolynomial_factors_destruct(sq_free_factors, 1);
  lp_upolynomial_delete(to_factor);

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_Zp("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(result, trace_out); tracef("\n");
  }

  return result;
}


/**
 * All polynomials here are in Zp for some prime p.
 *
 * Given a factorization A = A1*A2*...*Ar, with (Ai, Aj) = 1. Let Pi be the product
 * of all Aj but Ai, i.e. A/A_i, the function computes U1, ..., Ur such that
 *
 * [1]  P1*U1 + P2*U2 + ... + Pr*Ur = 1
 *
 * Note that above exists since gcd(P1, ..., Pn) = 1.
 *
 * The algorithm goes inductively. Let Qi = Ai*...*Ar and say we are solving as
 *
 * [2]  P1*U1 + P2*U2 + ... + Pr*Ur = D
 *
 * Then we can solve
 *
 * [3]  V1*A1 + U1*Q2 = D
 *
 * using extended gcd, obtaining U1, V1 with
 *
 * [4]  deg(V1) < deg(Q2) and deg(U1) < deg(A1).
 *
 * Replacing [3] in [2] we get (P1 = Q2)
 *
 * [5]  P2*U2 + ... + Pr*Ur = V1*A1
 *
 * Both sides are divisible by A1 producing an inductive subproblem
 *
 * [6] (P2/A1)*U2 + ... + (Pr/A1)*Ur = V1
 */
void hensel_lift_initialize(const lp_upolynomial_factors_t* A, lp_upolynomial_factors_t* U) {

  int i;

  if (trace_is_enabled("hensel")) {
    tracef("hensel_lift_initialize(");
    lp_upolynomial_factors_print(A, trace_out);
    tracef(")\n");
  }

  // The number of factors
  const int r = A->size;

  lp_int_ring_t* K = lp_upolynomial_factors_ring(A);

  // All the Q_i = A_i*...*A_r
  lp_upolynomial_t* Q[r];
  Q[r-1] = lp_upolynomial_construct_copy(A->factors[r-1]);
  for (i = r-2; i >= 0; -- i) {
    Q[i] = lp_upolynomial_mul(Q[i+1], A->factors[i]);
  }

  // The D we will be solving for
  lp_upolynomial_t* D = lp_upolynomial_construct_power(K, 0, 1);

  for (i = 0; i < r-1; ++ i) {
    lp_upolynomial_t* Ui = 0;
    lp_upolynomial_t* Vi = 0;

    // Solve D = Vi*Ai + U_i*Qi+1
    if (trace_is_enabled("hensel")) {
      tracef("A_i = "); lp_upolynomial_print(A->factors[i], trace_out); tracef("\n");
      tracef("Q_i+1 = "); lp_upolynomial_print(Q[i + 1], trace_out); tracef("\n");
      tracef("D  = "); lp_upolynomial_print(D, trace_out); tracef("\n");
    }

    lp_upolynomial_solve_bezout(A->factors[i], Q[i+1], D, &Vi, &Ui);

    if (trace_is_enabled("hensel")) {
      tracef("V_i = "); lp_upolynomial_print(Vi, trace_out); tracef("\n");
      tracef("U_i = "); lp_upolynomial_print(Ui, trace_out); tracef("\n");
    }

    // Remember the solution
    lp_upolynomial_factors_add(U, Ui, 1);

    // Next to solve for
    lp_upolynomial_delete(D);
    D = Vi;
  }

  // The last one
  lp_upolynomial_factors_add(U, D, 1);

  // Remove temporaries
  for (i = 0; i < r; ++ i) {
    lp_upolynomial_delete(Q[i]);
  }

  if (trace_is_enabled("hensel")) {
    tracef("hensel_lift_initialize(");
    lp_upolynomial_factors_print(A, trace_out);
    tracef(") => ");
    lp_upolynomial_factors_print(U, trace_out);
    tracef("\n");
  }
}

/**
 * Given the polynomials
 *
 *   A1, A2, ..., An in any K[x],
 *
 * the function computes the products
 *
 *   P1 =   *A2*A3*...*An
 *   P2 = A1*  *A3*...*An
 *   P3 = A1*A2*  *...*An
 *             ...
 *   Pn = A1*A2*A3*...*
 *
 * and returns the product
 *
 *   A1*A2*...*An.
 *
 * All products are in Z[x]
 */
void hensel_lift_compute_products(const lp_upolynomial_factors_t* A, lp_upolynomial_t** P_z) {

  int k;
  int n = A->size;

  assert(lp_upolynomial_factors_ring(A) != lp_Z);

  // Compute Z[x] lifts
  lp_upolynomial_t* A_z[n];
  for(k = 0; k < n; ++ k) {
    A_z[k] = lp_upolynomial_construct_copy_K(lp_Z, A->factors[k]);
  }

  // Compute P1[k] = A1*...*A{k-1}
  lp_upolynomial_t* P1_z[n];
  P1_z[0] = lp_upolynomial_construct_power(lp_Z, 0, 1);
  for(k = 1; k < n; ++ k) {
    P1_z[k] = lp_upolynomial_mul(P1_z[k-1], A_z[k-1]);
  }

  // Compute P2[k] = A{k+1}*...*An
  lp_upolynomial_t* P2_z[n];
  P2_z[n-1] = lp_upolynomial_construct_power(lp_Z, 0, 1);
  for(k = n-2; k >= 0; --k) {
    P2_z[k] = lp_upolynomial_mul(P2_z[k+1], A_z[k+1]);
  }

  // Compute the P[k]
  for(k = 0; k < n; ++ k) {
    P_z[k] = lp_upolynomial_mul(P1_z[k], P2_z[k]);
  }

  // Free the temps
  for(k = 0; k < n; ++ k) {
    lp_upolynomial_delete(A_z[k]);
    lp_upolynomial_delete(P1_z[k]);
    lp_upolynomial_delete(P2_z[k]);
  }

  if (trace_is_enabled("hensel")) {
    for (k = 0; k < n; ++ k) {
      tracef("P[%d] = ", k); lp_upolynomial_print(P_z[k], trace_out); tracef("\n");
    }
  }
}


/**
 * Perform a quadratic Hensel lift
 *
 * Input:
 *
 *  (1) q = p^k for some prime p
 *  (2) F: polynomial in Z[x]
 *  (3) A1, ..., Ar: polynomials in Zq[x]
 *  (4) U1, ..., Ur: polynomials in Zq[x]
 *
 * Assume:
 *
 *  [A1] F = A1 * A2 * ... * An            (mod q)
 *  [A2] gcd(lc(F), q) = 1,
 *       lc(F) = lc(A1) (mod q),
 *       lc(Ak) = 1 for k > 1,
 *
 * and with Pk = (A1*...*An)/Ak (product withouth Ak) assume
 *
 *  [A3] 1 = U1*P1 + U2*P2 + ... * Pn      (mod q)
 *  [A4] deg(Uk) < deg(Ak)
 *
 * Output:
 *
 *  Bk = Ak (mod q)
 *  Vk = Uk (mod q)
 *
 *  Same as input for (q^2, F, Bk, Vk)
 *
 */
void hensel_lift_quadratic(const lp_upolynomial_t* F,
    const lp_upolynomial_factors_t* A, const lp_upolynomial_factors_t* U,
    lp_upolynomial_factors_t* B, lp_upolynomial_factors_t* V)
{
  size_t k;

  assert(B->size == 0);
  assert(V->size == 0);

  // The ring we are lifting
  lp_int_ring_t* Zq = A->factors[0]->K;
  // The modulus
  const lp_integer_t* q = &Zq->M;

  // The lifted ring
  lp_integer_t qq;
  integer_construct_from_int(lp_Z, &qq, 0);
  integer_pow(lp_Z, &qq, q, 2);
  lp_int_ring_t* Zqq = lp_int_ring_create(&qq, 0);

  if (trace_is_enabled("hensel")) {
    tracef("hensel_lift_quadratic("); lp_upolynomial_print(F, trace_out); tracef(", ");
    lp_upolynomial_factors_print(A, trace_out);
    tracef(", "); lp_upolynomial_factors_print(U, trace_out); tracef(")\n");
    tracef("q   = "); integer_print(q, trace_out); tracef("\n");
    tracef("q^2 = "); integer_print(&qq, trace_out); tracef("\n");
  }

  /**
   * We are looking for
   *
   *   Bk = Ak (mod q)
   *   Bk = Ak + q*Sk
   *
   * such that
   *
   *   F = B1*B2*...*Br (mod q^2)
   *     = A1*A2*...*Ar + q*(S1*P1 + ... + Sr*Pr) (mod q^2)
   *
   * which is equivalent to
   *
   *   D = (F - A1*A2*...*Ar)/q = S1*P1 + ... + Sr*Pr (mod q)
   *
   * Here we know that deg(D) < deg(F) and therefore the equality holds mod
   * A1*A2*...*Ar also.
   *
   * Using [A3],[A4], we know that a general solutions are of the form
   *
   *   D*Uk
   *
   * but these do not satisfy degree requirements for the the lifting (deg(Sk) <
   * deg(Ak)). Since lc(Ak) is invertible in Zq, then
   *
   *   Sk = D*Uk (mod Ak) (in Zq[x])
   *
   * is also a solution.
   */

  lp_upolynomial_t* prod_Ak = lp_upolynomial_construct_copy_K(lp_Z, A->factors[0]);
  for (k = 1; k < A->size; ++ k) {
    lp_upolynomial_t* tmp1 = lp_upolynomial_construct_copy_K(lp_Z, A->factors[k]);
    lp_upolynomial_t* tmp2 = prod_Ak;
    prod_Ak = lp_upolynomial_mul(tmp1, tmp2);
    lp_upolynomial_delete(tmp1);
    lp_upolynomial_delete(tmp2);
  }

  if (trace_is_enabled("hensel")) {
    tracef("prod_Ak = "); lp_upolynomial_print(prod_Ak, trace_out); tracef("\n");
  }

  lp_upolynomial_t* F_sub_prod_Ak = lp_upolynomial_sub(F, prod_Ak);

  if (trace_is_enabled("hensel")) {
    tracef("F_sub_prod_Ak = "); lp_upolynomial_print(F_sub_prod_Ak, trace_out); tracef("\n");
  }

  lp_upolynomial_t* D = lp_upolynomial_div_exact_c(F_sub_prod_Ak, q);

  if (trace_is_enabled("hensel")) {
    tracef("D = "); lp_upolynomial_print(D, trace_out); tracef("\n");
  }

  lp_upolynomial_t* D_q = lp_upolynomial_construct_copy_K(Zq, D);

  if (trace_is_enabled("hensel")) {
    tracef("D_q = "); lp_upolynomial_print(D_q, trace_out); tracef("\n");
  }

  for (k = 0; k < A->size; ++ k) {

    lp_upolynomial_t* D_mult_Uk_q = lp_upolynomial_mul(D_q, U->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("D_mult_Uk_q = "); lp_upolynomial_print(D_mult_Uk_q, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Sk_q = lp_upolynomial_rem_exact(D_mult_Uk_q, A->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("Sk_q = "); lp_upolynomial_print(Sk_q, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Sk_qq = lp_upolynomial_construct_copy_K(Zqq, Sk_q);
    if (trace_is_enabled("hensel")) {
      tracef("Sk_qq = "); lp_upolynomial_print(Sk_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Ak_qq = lp_upolynomial_construct_copy_K(Zqq, A->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("Ak_qq = "); lp_upolynomial_print(Ak_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Sk_times_q_qq = lp_upolynomial_mul_c(Sk_qq, q);
    if (trace_is_enabled("hensel")) {
      tracef("Sk_times_q_qq = "); lp_upolynomial_print(Sk_times_q_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Bk_qq = lp_upolynomial_add(Ak_qq, Sk_times_q_qq);
    if (trace_is_enabled("hensel")) {
      tracef("Bk_qq = "); lp_upolynomial_print(Bk_qq, trace_out); tracef("\n");
    }

    lp_upolynomial_factors_add(B, Bk_qq, 1);

    lp_upolynomial_delete(D_mult_Uk_q);
    lp_upolynomial_delete(Sk_q);
    lp_upolynomial_delete(Sk_qq);
    lp_upolynomial_delete(Ak_qq);
    lp_upolynomial_delete(Sk_times_q_qq);
  }

  /**
   * In the same manner as above we are looking for
   *
   *  Vk = Uk (mod q)
   *  Vk = Uk + q*Tk
   *
   * such that
   *
   *  1 = V1*P1 + V2*P2 + ... + Vr*Pr (mod q^2)
   *    = (U1 + q*T1)*P1 + ...        (mod q^2)
   *
   * which is equivalent to
   *
   *  E = (1 - (U1*P1 + ... + Ur*Pr))/q = T1*P1 + ... + Tr*Pr
   *
   * again, with deg(E) < deg(F) so the equality also holds mod A1*...*Ar.
   *
   * Using [A3] again, the solution is of the form
   *
   *  E*Uk
   *
   * but we don't satisfy degree requirements. Again, we compute modulo Ak
   * obtaining
   *
   *  Tk = E*Uk (mod Ak) (in Zq[x])
   */
  lp_upolynomial_t* P[A->size];
  hensel_lift_compute_products(B, P);
  lp_upolynomial_t* E = lp_upolynomial_construct_power(lp_Z, 0, 1);
  if (trace_is_enabled("hensel")) {
    tracef("E = "); lp_upolynomial_print(E, trace_out); tracef("\n");
  }
  for (k = 0; k < A->size; ++ k) {
    lp_upolynomial_t* Uk = lp_upolynomial_construct_copy_K(lp_Z, U->factors[k]);
    lp_upolynomial_t* mul = lp_upolynomial_mul(Uk, P[k]);
    lp_upolynomial_t* sub = lp_upolynomial_sub(E, mul);
    lp_upolynomial_delete(Uk);
    lp_upolynomial_delete(mul);
    lp_upolynomial_delete(E);
    E = sub;
    if (trace_is_enabled("hensel")) {
      tracef("E = "); lp_upolynomial_print(E, trace_out); tracef("\n");
    }
  }
  lp_upolynomial_t* tmp = E;
  E = lp_upolynomial_div_exact_c(E, q);
  if (trace_is_enabled("hensel")) {
    tracef("E = "); lp_upolynomial_print(E, trace_out); tracef("\n");
  }
  lp_upolynomial_delete(tmp);
  lp_upolynomial_t* E_q = lp_upolynomial_construct_copy_K(Zq, E);
  if (trace_is_enabled("hensel")) {
    tracef("E_q = "); lp_upolynomial_print(E_q, trace_out); tracef("\n");
  }

  for (k = 0; k < A->size; ++ k) {

    lp_upolynomial_t* E_times_Uk_q = lp_upolynomial_mul(E_q, U->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("E_times_Uk_q = "); lp_upolynomial_print(E_times_Uk_q, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Tk_q = lp_upolynomial_rem_exact(E_times_Uk_q, A->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("Tk_q = "); lp_upolynomial_print(Tk_q, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Tk_qq = lp_upolynomial_construct_copy_K(Zqq, Tk_q);
    if (trace_is_enabled("hensel")) {
      tracef("Tk_qq = "); lp_upolynomial_print(Tk_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Tk_times_q_qq = lp_upolynomial_mul_c(Tk_qq, q);
    if (trace_is_enabled("hensel")) {
      tracef("Tk_times_q_qq = "); lp_upolynomial_print(Tk_times_q_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Uk_qq = lp_upolynomial_construct_copy_K(Zqq, U->factors[k]);
    if (trace_is_enabled("hensel")) {
      tracef("Uk_qq = "); lp_upolynomial_print(Uk_qq, trace_out); tracef("\n");
    }
    lp_upolynomial_t* Vk_qq = lp_upolynomial_add(Uk_qq, Tk_times_q_qq);
    if (trace_is_enabled("hensel")) {
      tracef("Vk_qq = "); lp_upolynomial_print(Vk_qq, trace_out); tracef("\n");
    }

    lp_upolynomial_factors_add(V, Vk_qq, 1);

    lp_upolynomial_delete(E_times_Uk_q);
    lp_upolynomial_delete(Tk_q);
    lp_upolynomial_delete(Tk_qq);
    lp_upolynomial_delete(Tk_times_q_qq);
    lp_upolynomial_delete(Uk_qq);
  }

  if (trace_is_enabled("hensel")) {
    tracef("Lifted to:");
    tracef("\nB = "); lp_upolynomial_factors_print(B, trace_out);
    tracef("\nV =  "); lp_upolynomial_factors_print(V, trace_out);
    tracef("\n");
  }

  // Remove the temporaries
  integer_destruct(&qq);
  lp_int_ring_detach(Zqq);
  for (k = 0; k < A->size; ++ k) {
    lp_upolynomial_delete(P[k]);
  }
  lp_upolynomial_delete(prod_Ak);
  lp_upolynomial_delete(F_sub_prod_Ak);
  lp_upolynomial_delete(D);
  lp_upolynomial_delete(D_q);
  lp_upolynomial_delete(E);
  lp_upolynomial_delete(E_q);
}



/**
 * Reconstruct factorization of f, from it's lifted factorization factors_p. Put
 * the result in factors. Basically try combinations of factors and see if they
 * divide f.
 */
void factorization_recombination(const lp_upolynomial_t* f, const lp_upolynomial_factors_t* factors_p, lp_upolynomial_factors_t* factors) {

  // Our copy of f to work with
  lp_upolynomial_t* to_factor = lp_upolynomial_construct_copy(f);

  // Some counters
  int i, k;

  // Maximal degree we want to consider
  size_t max_degree = lp_upolynomial_degree(f)/2;

  const int size = factors_p->size;

  // Factors we are still considering (initially all with degree <= max_degree)
  int enabled_count = 0;
  int enabled[size+1];
  for (k = 0; k < size; ++ k) {
    if (lp_upolynomial_degree(factors_p->factors[k]) <= max_degree) {
      enabled[k] = 1;
      enabled_count ++;
    } else {
      enabled[k] = 0;
    }
  }
  enabled[size] = 1; // max is always enabled

  // Current selection that we will be updating (-1 for not used)
  int sel[size+1];
  for (k = 0; k < size+1; ++ k) {
    sel[k] = size; // All should be smaller than size
  }

  // Each loop iteration finds a new combination. At each iteration we know that
  // we have exhausted all combinations of the same length that are smaller or
  // equal to the current one.
  //
  // At each step, if a selection is in a non-enabled state, it was just
  // eliminated. Therefore the whole selection needs to be updated.
  int done = 0;
  int sel_size = 1;
  while (!done) {

    int carry = 1;   // Are we carrying forward
    int minimal = 0; // We are minimal up to (not including) this position

    for (k = 0; (carry || !enabled[sel[k]]); ++ k) {

      // If out of bound, we're done
      if (k == size) {
        done = 1;
        break;
      }

      // If we are increasing the next one, the preceding sequence
      // needs to be set to minimal
      for (i = minimal; i < k; ++ i) {
        sel[i] = i > 0 ? sel[i-1] + 1 : 0;
        while (!enabled[sel[i]]) {
          sel[i]++;
        }
        // If some of the selection is the max, we're done
        if (sel[i] == size) {
          done = 1;
          break;
        }
        // We are now minimal up to (not including) k
        minimal = k;
      }

      // Are we done
      if (done) {
        break;
      }

      // Initialize the current one to the first possible value
      // * not the mark for max
      // * at least the previous one
      if (k > 0) {
        if (sel[k] == size || sel[k] < sel[k-1]) {
          sel[k] = sel[k-1] + 1;
          carry = 0;
        }
      } else if (sel[k] == size) {
        sel[k] = 0;
        carry = 0;
      }

      // Increase until ok
      while (sel[k] < sel[k+1] && (carry || !enabled[sel[k]])) {
        sel[k] ++;
        carry = 0;
      }

      // Reached top, so we increase the size
      if (sel[k] == size) {
        carry = 1;
        sel_size ++;
      }

      // If we made it, no need for more carry
      carry = !(sel[k] < sel[k+1]);
    }

    // If we are not done, try this combination of factors
    if (!done) {

      if (!done && trace_is_enabled("factorization")) {
        tracef("sel = [");
        for (i = 0; i < sel_size; ++ i) {
          if (i) tracef(" ");
          tracef("%d", sel[i]);
        }
        tracef("]\n");
      }

      // Check the sum of degrees
      size_t deg_sum = 0;
      for (i = 0; i < sel_size; ++ i) {
        deg_sum += lp_upolynomial_degree(factors_p->factors[sel[i]]);
      }

      if (deg_sum <= max_degree) {

        // Construct the candidate
        lp_upolynomial_t* candidate = lp_upolynomial_construct_copy(factors_p->factors[sel[0]]);
        for (i = 1; i < sel_size; ++ i) {
          lp_upolynomial_t* tmp = candidate;
          candidate = lp_upolynomial_mul(candidate, factors_p->factors[sel[i]]);
          lp_upolynomial_delete(tmp);
        }
        lp_upolynomial_set_ring(candidate, lp_Z);
        if (trace_is_enabled("factorization")) {
          tracef("candidate = "); lp_upolynomial_print(candidate, trace_out); tracef("\n");
        }

        // Check divisibility
        if (lp_upolynomial_divides(candidate, to_factor)) {
          // Hooray, we found a factor
          lp_upolynomial_factors_add(factors, candidate, 1);
          lp_upolynomial_t* tmp = to_factor;
          to_factor = lp_upolynomial_div_exact(to_factor, candidate);
          lp_upolynomial_delete(tmp);
          // Disable current selection
          for (i = 0; i < sel_size; ++ i) {
            enabled[sel[i]] = 0;
          }
        }
      }
    }
  }

  // Last one standing
  if (lp_upolynomial_degree(to_factor) == 0) {
    lp_upolynomial_delete(to_factor);
  } else {
    lp_upolynomial_factors_add(factors, to_factor, 1);
  }
}

// First 100 primes
static const long primes[100] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
   31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
   73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
  283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
  353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
  419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
  467, 479, 487, 491, 499, 503, 509, 521, 523, 541
};

static size_t primes_count = sizeof(primes)/sizeof(long);

/**
 * Factors a square-free f in Z.
 */
lp_upolynomial_factors_t* upolynomial_factor_Z_square_free(const lp_upolynomial_t* f) {

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_Z_square_free("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }

  assert(f->K == lp_Z);
  assert(lp_upolynomial_degree(f) > 1);

  // The result
  lp_upolynomial_factors_t* factors = lp_upolynomial_factors_construct();

  // Get the bound on coefficients sufficient when lifting to Z
  lp_integer_t coefficient_bound;
  integer_construct_from_int(lp_Z, &coefficient_bound, 0);
  upolynomial_factor_bound_landau_mignotte(f, lp_upolynomial_degree(f)/2, &coefficient_bound);
  integer_mul_int(lp_Z, &coefficient_bound, &coefficient_bound, 2);

  if (trace_is_enabled("factorization")) {
    tracef("coefficient_bound = "); integer_print(&coefficient_bound, trace_out); tracef("\n");
  }

  // The prime factorization
  lp_upolynomial_factors_t* factors_p_best = 0;

  // Find a good prime modulus
  size_t prime_i;
  int candidates_to_check = 1;
  for (prime_i = 0; candidates_to_check && prime_i < primes_count; ++ prime_i) {

    // The prime we'll try
    lp_integer_t prime;
    integer_construct_from_int(lp_Z, &prime, primes[prime_i]);

    if (trace_is_enabled("factorization")) {
      tracef("prime = "); integer_print(&prime, trace_out); tracef("\n");
    }

    // We need a prime that keeps the degree of f, and maintains being square-free
    if (!integer_divides(lp_Z, &prime, lp_upolynomial_lead_coeff(f))) {

      // The ring to use
      lp_int_ring_t* K_p = lp_int_ring_create(&prime, 1);

      // compute GCD of (f, f') in Z_p (to check for square-free
      lp_upolynomial_t* f_p = lp_upolynomial_construct_copy_K(K_p, f);
      if (trace_is_enabled("factorization")) {
        tracef("f_p = "); lp_upolynomial_print(f_p, trace_out); tracef("\n");
      }
      lp_upolynomial_t* f_d_p = lp_upolynomial_derivative(f_p);
      lp_upolynomial_t* gcd_p = lp_upolynomial_gcd(f_p, f_d_p);

      if (lp_upolynomial_is_one(gcd_p)) {

        // Initial factorization in in Z_p
        // We know that factorization will be square-free due to gcd = 1
        lp_upolynomial_factors_t* factors_p = lp_upolynomial_factor(f_p);

        // If there is more than one factor we do the lifting
        if (!factors_p_best || factors_p->size < factors_p_best->size) {
          if (factors_p_best) {
            lp_upolynomial_factors_destruct(factors_p_best, 1);
          }
          factors_p_best = factors_p;
        }

        if (factors_p_best->size == 1) {
          // We're done
          candidates_to_check = 0;
        } else {
          // One less candidate
          candidates_to_check --;
        }
      }

      // Get rid of temps
      lp_upolynomial_delete(f_p);
      lp_upolynomial_delete(f_d_p);
      lp_upolynomial_delete(gcd_p);
      lp_int_ring_detach(K_p);
    }

    integer_destruct(&prime);
  }

  assert(factors_p_best);

  if (factors_p_best->size > 1) {
    // The U's for the Hensel lift
    lp_upolynomial_factors_t* U_p = lp_upolynomial_factors_construct();
    hensel_lift_initialize(factors_p_best, U_p);
    lp_upolynomial_factors_t* factors_q = lp_upolynomial_factors_construct();
    lp_upolynomial_factors_t* U_q = lp_upolynomial_factors_construct();

    // Lift the factorization until enough to reconstruct in Z
    while (integer_cmp(lp_Z, &lp_upolynomial_factors_ring(factors_p_best)->M,
        &coefficient_bound) < 0) {
      // Do the lift
      hensel_lift_quadratic(f, factors_p_best, U_p, factors_q, U_q);
      // Swap and clear
      lp_upolynomial_factors_swap(factors_p_best, factors_q);
      lp_upolynomial_factors_swap(U_p, U_q);
      lp_upolynomial_factors_clear(factors_q);
      lp_upolynomial_factors_clear(U_q);
    }

    // Do the reconstruction
    factorization_recombination(f, factors_p_best, factors);

    // Remove temps
    lp_upolynomial_factors_destruct(U_p, 1);
    lp_upolynomial_factors_destruct(U_q, 1);
    lp_upolynomial_factors_destruct(factors_q, 1);
  } else {
    // Primitive
    lp_upolynomial_factors_add(factors, lp_upolynomial_construct_copy(f), 1);
  }

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_Z_square_free("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  // Free the temps
  lp_upolynomial_factors_destruct(factors_p_best, 1);
  integer_destruct(&coefficient_bound);

  // Done
  return factors;
}

lp_upolynomial_factors_t* upolynomial_factor_Z(const lp_upolynomial_t* f) {

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_Z("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }

  assert(f->K == lp_Z);

  // The result
  lp_upolynomial_factors_t* factors = lp_upolynomial_factors_construct();

  // The content of the polynomial (constant of the factorization)
  lp_upolynomial_content_Z(f, &factors->constant);

  // We factor the primitive part
  lp_upolynomial_t* f_pp = lp_upolynomial_primitive_part_Z(f);

  // Get a square-free decomposition of f
  lp_upolynomial_factors_t* sq_free_factors = lp_upolynomial_factor_square_free_primitive(f_pp);
  assert(integer_cmp_int(lp_Z, &sq_free_factors->constant, 1) == 0);

  // Factor individuals
  size_t sq_free_factor_i;
  for (sq_free_factor_i = 0; sq_free_factor_i < sq_free_factors->size; ++ sq_free_factor_i) {
    lp_upolynomial_t* sq_free_factor = sq_free_factors->factors[sq_free_factor_i];
    size_t sq_free_factor_multiplicity = sq_free_factors->multiplicities[sq_free_factor_i];
    size_t sq_free_factor_deg = lp_upolynomial_degree(sq_free_factor);

    // Do special cases
    if (sq_free_factor_deg == 1) {
      lp_upolynomial_factors_add(factors, sq_free_factor, sq_free_factor_multiplicity);
      continue;
    }

    // Factorize the square-free factor
    lp_upolynomial_factors_t* sq_free_factor_factors = upolynomial_factor_Z_square_free(sq_free_factor);

    // We don't need it anymore
    lp_upolynomial_delete(sq_free_factor);

    // Add the factors
    size_t current;
    for (current = 0; current < sq_free_factor_factors->size; ++ current) {
      lp_upolynomial_t* factor = sq_free_factor_factors->factors[current];
      size_t factor_multiplicity = sq_free_factor_factors->multiplicities[current];
      lp_upolynomial_factors_add(factors, factor, sq_free_factor_multiplicity * factor_multiplicity);
    }

    // Remove the factorization
    lp_upolynomial_factors_destruct(sq_free_factor_factors, 0);
  }

  // Get rid of the square free factors
  lp_upolynomial_factors_destruct(sq_free_factors, 0);

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_Z("); lp_upolynomial_print(f, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  // Done
  return factors;
}
