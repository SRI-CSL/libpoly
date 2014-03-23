/*
 * univariate_polynomial_gcd.c
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#include "upolynomial/gcd.h"
#include "upolynomial/upolynomial_dense.h"
#include "upolynomial/upolynomial.h"

#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <assert.h>

STAT_DECLARE(int, upolynomial, gcd_euclid);
STAT_DECLARE(int, upolynomial, gcd_euclid_extended);
STAT_DECLARE(int, upolynomial, gcd_subresultant);
STAT_DECLARE(int, upolynomial, gcd_heuristic);
STAT_DECLARE(int, upolynomial, gcd_heuristic_success);

/**
 * Computing using Euclid's algorithm.
 *
 * For regular GCD we have polynomials A, B, we set
 *
 *   r_0 = A,
 *   r_1 = B,
 *
 * and we go on by computing
 *
 *   r_0 = r_1*q + r2
 *
 * then
 *
 *   gcd(A, B) = gcd(r_0, r_1) = gcd(r_1, r_2)
 *
 * In addition if we take out the content of A and B out in the beginning, we
 * are also sure that
 *
 *   gcd(A, B) = gcd(r_1, pp(r_2)).
 *
 * To compute the extended GCD we start with
 *
 *   r_0 = A    r_1 = B
 *   s_0 = 1    s_1 = 0
 *   t_0 = 0    t_1 = 1
 *
 * We maintain the invariants that
 *
 * [I]  r_i = A*s_i + B*t_i.
 *
 * We are done r_n = 0, r_{n-1} = gcd, so we take s_{n-1}, t_{n-1} as solutions.
 *
 *   q   = r_0 / r_1
 *   r_2 = r_0 - q*r_1
 *   s_2 = s_0 - q*s_1
 *   t_2 = t_0 - q*t_1
 *
 * Each pair gets multiplied by (1, -q) so it's easy to show that the
 * invariant [I] holds.
 */
upolynomial_t* upolynomial_gcd_euclid(const upolynomial_t* A, const upolynomial_t* B, upolynomial_t** U, upolynomial_t** V)
{
  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("upolynomial_gcd_euclid(%P, %P)\n", A, B);
  }

  assert(!upolynomial_ops.is_zero(B));

  int extended_gcd = (U != 0 && V != 0);

  if (extended_gcd) STAT(upolynomial, gcd_euclid_extended) ++;
  else STAT(upolynomial, gcd_euclid) ++;

  // Degrees of p and q
  size_t deg_A = upolynomial_ops.degree(A);
  assert(deg_A >= upolynomial_ops.degree(B));

  // The ring of computation
  assert(A->K == B->K);
  int_ring K = A->K;
  assert(K && K->is_prime);

  // The remainder in the calculation
  upolynomial_t* D = 0;

  // Buffers to keep A, B, and remainders
  upolynomial_dense_t r_0;
  upolynomial_dense_t r_1;
  upolynomial_dense_ops.construct_p(&r_0, deg_A + 1, A);
  upolynomial_dense_ops.construct_p(&r_1, deg_A + 1, B);

  // Buffers to keep div and rem
  upolynomial_dense_t r_2;
  upolynomial_dense_t q;
  upolynomial_dense_ops.construct(&r_2, deg_A + 1); // Some extra for reuse
  upolynomial_dense_ops.construct(&q, deg_A + 1);

  // Buffers for the extended gcd computation
  upolynomial_dense_t s_0, s_1;
  upolynomial_dense_t t_0, t_1;
  if (extended_gcd) {
    upolynomial_dense_ops.construct(&s_0, deg_A + 1);
    upolynomial_dense_ops.construct(&s_1, deg_A + 1);
    upolynomial_dense_ops.construct(&t_0, deg_A + 1);
    upolynomial_dense_ops.construct(&t_1, deg_A + 1);
    // s_0, t_1 = 1
    integer_ops.assign_int(Z, s_0.coefficients, 1);
    integer_ops.assign_int(Z, t_1.coefficients, 1);
  }

  do {
    // One step of division
    upolynomial_dense_ops.div_general(K, 1 /* exact */, &r_0, &r_1, &q, &r_2);

    if (debug_trace_ops.is_enabled("gcd")) {
      tracef("r_0 = ");
      upolynomial_dense_ops.print(&r_0, trace_out);
      tracef("\nr_1 = ");
      upolynomial_dense_ops.print(&r_1, trace_out);
      tracef("\nq = ");
      upolynomial_dense_ops.print(&q, trace_out);
      tracef("\nr_2 = ");
      upolynomial_dense_ops.print(&r_2, trace_out);
      if (extended_gcd) {
        tracef("\ns_0 = ");
        upolynomial_dense_ops.print(&s_0, trace_out);
        tracef("\ns_1 = ");
        upolynomial_dense_ops.print(&s_1, trace_out);
        tracef("\nt_0 = ");
        upolynomial_dense_ops.print(&t_0, trace_out);
        tracef("\nt_1 = ");
        upolynomial_dense_ops.print(&t_1, trace_out);
      }
      tracef("\n");
    }

    // Check we are done
    if (upolynomial_dense_ops.is_zero(&r_2))  {
      // We're in a field, make it monic
      integer_t lc;
      integer_ops.construct_copy(K, &lc, r_1.coefficients + r_1.size - 1);
      if (integer_ops.cmp_int(Z, &lc, 1)) {
        upolynomial_dense_ops.div_c(&r_1, K, &lc);
        if (extended_gcd) {
          upolynomial_dense_ops.div_c(&s_1, K, &lc);
          upolynomial_dense_ops.div_c(&t_1, K, &lc);
        }
      }
      integer_ops.destruct(&lc);
      D = upolynomial_dense_ops.to_upolynomial(&r_1, K);
      if (extended_gcd) {
        *U = upolynomial_dense_ops.to_upolynomial(&s_1, K);
        *V = upolynomial_dense_ops.to_upolynomial(&t_1, K);
      }
    } else {

      // Extended gcd computation
      if (extended_gcd) {
        // s2 = s0 - q*s1
        // (s0, s1) = (s1, s2)
        upolynomial_dense_ops.sub_mult(&s_0, K, &q, &s_1);
        upolynomial_dense_ops.swap(&s_0, &s_1);;
        // t2 = t0 - q*t2
        upolynomial_dense_ops.sub_mult(&t_0, K, &q, &t_1);
        upolynomial_dense_ops.swap(&t_0, &t_1);
      }

      // (r0, r1) = (r1, r2)
      upolynomial_dense_ops.swap(&r_0, &r_1);
      upolynomial_dense_ops.swap(&r_1, &r_2);


    }

  } while (D == 0);

  if (extended_gcd) {
    upolynomial_dense_ops.destruct(&s_0);
    upolynomial_dense_ops.destruct(&s_1);
    upolynomial_dense_ops.destruct(&t_0);
    upolynomial_dense_ops.destruct(&t_1);
  }

  upolynomial_dense_ops.destruct(&q);
  upolynomial_dense_ops.destruct(&r_2);
  upolynomial_dense_ops.destruct(&r_0);
  upolynomial_dense_ops.destruct(&r_1);

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("upolynomial_gcd_euclid(%P, %P) = %P\n", A, B, D);
  }

  return D;
}

upolynomial_t* upolynomial_gcd_subresultant(const upolynomial_t* A, const upolynomial_t* B) {

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("upolynomial_gcd_subresultant(%P, %P)\n", A, B);
  }
  STAT(upolynomial, gcd_subresultant) ++;

  assert(!upolynomial_ops.is_zero(B));

  // The ring of compuation
  assert(A->K == B->K);
  assert(A->K == Z);
  int_ring K = A->K;

  // Degrees of p and q
  size_t deg_A = upolynomial_ops.degree(A);
  assert(deg_A >= upolynomial_ops.degree(B));

  // The remainder in the calculation
  upolynomial_t* D = 0;

  // Dense representations of the remainders to keep p and q
  upolynomial_dense_t r_0;
  upolynomial_dense_t r_1;
  upolynomial_dense_ops.construct_p(&r_0, deg_A + 1, A);
  upolynomial_dense_ops.construct_p(&r_1, deg_A + 1, B);

  // Contents of a, b
  integer_t A_cont, B_cont;
  integer_ops.construct_from_int(K, &A_cont, 0);
  integer_ops.construct_from_int(K, &B_cont, 0);
  upolynomial_ops.content_Z(A, &A_cont);
  upolynomial_ops.content_Z(B, &B_cont);

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("cont(p) = "); integer_print(&A_cont, trace_out); tracef("\n");
    tracef("cont(q) = "); integer_print(&B_cont, trace_out); tracef("\n");
  }

  // d = gcd(content(p), content(q)))
  integer_t d;
  integer_ops.construct_from_int(K, &d, 1);

  integer_ops.gcd_Z(&d, &A_cont, &B_cont);
  if (integer_ops.cmp_int(Z, &d, 1)) {
    // GCD != 1
    upolynomial_dense_ops.div_c(&r_0, K, &A_cont);
    upolynomial_dense_ops.div_c(&r_1, K, &B_cont);
  }

  // Buffers to keep div and rem
  upolynomial_dense_t r_2;
  upolynomial_dense_t q;
  upolynomial_dense_ops.construct(&r_2, deg_A + 1); // Some extra for reuse
  upolynomial_dense_ops.construct(&q, deg_A + 1);

  // Adjustment coefficients
  integer_t g, h;
  integer_ops.construct_from_int(Z, &g, 1);
  integer_ops.construct_from_int(Z, &h, 1);

  // Temps for computation
  integer_t tmp1, tmp2;
  integer_ops.construct_from_int(K, &tmp1, 0);
  integer_ops.construct_from_int(K, &tmp2, 0);

  do {
    // \delta = deg(p) - deg(q)
    int delta = r_0.size - r_1.size;

    // One step of division
    upolynomial_dense_ops.div_general(K, 0, &r_0, &r_1, &q, &r_2);

    if (debug_trace_ops.is_enabled("gcd")) {
      tracef("r_0 = ");
      upolynomial_dense_ops.print(&r_0, trace_out);
      tracef("\nr_q = ");
      upolynomial_dense_ops.print(&r_1, trace_out);
      tracef("\nq = ");
      upolynomial_dense_ops.print(&q, trace_out);
      tracef("\nr_w = ");
      upolynomial_dense_ops.print(&r_2, trace_out);
      tracef("\n");
    }

    // Check if the remainder is of degree 0
    if (r_2.size == 1)  {
      if (integer_ops.sgn(Z, r_2.coefficients)) {
        // rem != 0, GCD(p, q) is 1 => total gcd is d
        integer_ops.assign(K, r_2.coefficients, &d);
        D = upolynomial_dense_ops.to_upolynomial(&r_2, K);
      } else {
        // rem == 0
        upolynomial_dense_ops.mk_primitive_Z(&r_1, 1);
        upolynomial_dense_ops.mult_c(&r_1, K, &d);
        D = upolynomial_dense_ops.to_upolynomial(&r_1, K);
      }
    } else {
      // p = q
      upolynomial_dense_ops.swap(&r_0, &r_1);
      // q = rem/(g*(h^delta)
      integer_ops.pow(K, &tmp1, &h, delta);
      integer_ops.mul(K, &tmp2, &tmp1, &g);
      upolynomial_dense_ops.swap(&r_1, &r_2);
      upolynomial_dense_ops.div_c(&r_1, K, &tmp2);
      // g = lc(p)
      integer_ops.assign(K, &g, r_0.coefficients + r_0.size - 1);
      // h = h^(1-delta)*g^delta = g^delta/(h^(delta-1))
      integer_ops.pow(K, &tmp1, &g, delta);
      integer_ops.pow(K, &tmp2, &h, delta - 1);
      integer_ops.div_exact(K, &h, &tmp1, &tmp2);
    }
  } while (D == 0);

  integer_ops.destruct(&tmp1);
  integer_ops.destruct(&tmp2);
  integer_ops.destruct(&g);
  integer_ops.destruct(&h);

  integer_ops.destruct(&d);
  integer_ops.destruct(&A_cont);
  integer_ops.destruct(&B_cont);

  upolynomial_dense_ops.destruct(&q);
  upolynomial_dense_ops.destruct(&r_2);
  upolynomial_dense_ops.destruct(&r_0);
  upolynomial_dense_ops.destruct(&r_1);

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("upolynomial_gcd_subresultant(%P, %P) = %P\n", A, B, D);
  }

  return D;
}

static void evaluate_polynomial(const upolynomial_t* A, const integer_t* A_content, unsigned pow, integer_t* out) {

  assert(pow > 0);

  integer_ops.assign_int(Z, out, 0);

  integer_t add;
  integer_t coeff;
  integer_ops.construct_from_int(Z, &add, 0);
  integer_ops.construct_from_int(Z, &coeff, 0);

  int k;
  for (k = 0; k < A->size; ++ k) {
    integer_ops.div_exact(Z, &coeff, &A->monomials[k].coefficient, A_content);
    if (A->monomials[k].degree == 0) {
      integer_ops.assign(Z, &add, &coeff);
    } else {
      integer_ops.mul_pow2(Z, &add, &coeff, A->monomials[k].degree * pow);
    }
    integer_ops.add(Z, out, out, &add);

    if (debug_trace_ops.is_enabled("gcd")) {
      tracef("out = "); integer_print(out, trace_out); tracef("\n");
    }
  }

  integer_ops.destruct(&coeff);
  integer_ops.destruct(&add);
}

/**
 * Reconstruct the polynomial from the value at 2^power, while multuplying with cont.
 */
static upolynomial_t* reconstruct_polynomial(size_t max_size, integer_t* p_value, unsigned pow, const integer_t* cont) {

  assert(pow > 0);

  // computation here
  upolynomial_dense_t p_d;
  upolynomial_dense_ops.construct(&p_d, max_size);

  integer_t div;
  integer_t rem;
  integer_ops.construct_from_int(Z, &div, 0);
  integer_ops.construct_from_int(Z, &rem, 0);

  integer_t P;
  integer_ops.construct_from_int(Z, &P, 2);
  integer_ops.pow(Z, &P, &P, pow);

  // Basically compute digits modulo x_value
  int d = 0;
  while (integer_ops.sgn(Z, p_value)) {
    integer_ops.div_rem_pow2_Z(&div, &rem, p_value, pow);

    if (integer_ops.bits(&rem) + 1 >= pow) {
      // If biger than 2^n/2 then take subtract 2^n-1
      integer_ops.sub(Z, &rem, &rem, &P);
      integer_ops.inc(Z, &div);
    }

    integer_ops.swap(Z, &rem, p_d.coefficients + d);
    integer_ops.swap(Z, &div, p_value);

    d ++;
  }
  upolynomial_dense_ops.touch(&p_d, Z, d - 1);

  // We only care about primitive GCDs in the reconstruction
  upolynomial_dense_ops.mk_primitive_Z(&p_d, 1);

  // Now, multiply with cont
  upolynomial_dense_ops.mult_c(&p_d, Z, cont);

  // Get the sparse representation
  upolynomial_t* p = upolynomial_dense_ops.to_upolynomial(&p_d, Z);

  // Free temporaries
  upolynomial_dense_ops.destruct(&p_d);
  integer_ops.destruct(&div);
  integer_ops.destruct(&rem);
  integer_ops.destruct(&P);

  return p;
}

int bound_valuation(const upolynomial_t* A, const upolynomial_t* B, const integer_t* A_cont, const integer_t* B_cont) {
  // Get the highest power 2^n such that
  //  2^n >= 2*min(|pp(A)|, |pp(B)|) + 2

  int A_max = 0;
  int B_max = 0;
  int k;

  integer_t tmp;
  integer_ops.construct_from_int(Z, &tmp, 0);

  for (k = 0; k < A->size; ++ k) {
    integer_ops.div_Z(&tmp, &A->monomials[k].coefficient, A_cont);
    int bits = integer_ops.bits(&tmp);
    if (bits > A_max) {
      A_max = bits;
    }
  }

  for (k = 0; k < B->size; ++ k) {
    integer_ops.div_Z(&tmp, &B->monomials[k].coefficient, B_cont);
    int bits = integer_ops.bits(&tmp);
    if (bits > B_max) {
      B_max = bits;
    }
  }

  integer_ops.destruct(&tmp);

  int bits_min = A_max > B_max ? B_max : A_max;

  return bits_min + 2;
}

upolynomial_t* upolynomial_gcd_heuristic(const upolynomial_t* A, const upolynomial_t* B, int attempts) {

  // Let's keep the smaller one in B
  if (upolynomial_ops.degree(A) < upolynomial_ops.degree(B)) {
    return upolynomial_gcd_heuristic(B, A, attempts);
  }

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("upolynomial_gcd_heuristic(%P, %P)\n", A, B);
  }
  STAT(upolynomial, gcd_heuristic) ++;

  upolynomial_t* D = 0;

  // The ring of computation
  assert(A->K == B->K);
  assert(A->K == Z);

  // content(A), content(B)
  integer_t A_cont, B_cont;
  integer_ops.construct_from_int(Z, &A_cont, 0);
  integer_ops.construct_from_int(Z, &B_cont, 0);
  upolynomial_ops.content_Z(A, &A_cont);
  upolynomial_ops.content_Z(B, &B_cont);

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("cont(p) = "); integer_print(&A_cont, trace_out); tracef("\n");
    tracef("cont(q) = "); integer_print(&B_cont, trace_out); tracef("\n");
  }

  // d = gcd(content(A), content(B)))
  integer_t d;
  integer_ops.construct_from_int(Z, &d, 1);
  integer_ops.gcd_Z(&d, &A_cont, &B_cont);

  if (debug_trace_ops.is_enabled("gcd")) {
    tracef("d = "); integer_print(&d, trace_out); tracef("\n");
  }

  // The number we use for valuation 2^n
  int n = bound_valuation(A, B, &A_cont, &B_cont);

  integer_t A_v, B_v, D_v;
  integer_ops.construct_from_int(Z, &A_v, 0);
  integer_ops.construct_from_int(Z, &B_v, 0);
  integer_ops.construct_from_int(Z, &D_v, 0);

  while (D == 0 && (attempts --)) {

    // Evaluate A and B
    evaluate_polynomial(A, &A_cont, n, &A_v);
    evaluate_polynomial(B, &B_cont, n, &B_v);

    if (debug_trace_ops.is_enabled("gcd")) {
        tracef("value of A/cont = "); integer_print(&A_v, trace_out); tracef("\n");
        tracef("value of B/cont = "); integer_print(&B_v, trace_out); tracef("\n");
    }

    // Get the gcd of the values and reconstruct the possible gcd
    integer_ops.gcd_Z(&D_v, &A_v, &B_v);
    // This also changes the value D_v
    D = reconstruct_polynomial(upolynomial_ops.degree(A) + 1, &D_v, n, &d);

    // Check if it divides both (B is smaller degree)
    if (!upolynomial_ops.divides(D, B) || !upolynomial_ops.divides(D, A)) {
      upolynomial_ops.destruct(D);
      D = 0;
    }

    // Try next power
    n ++;
  }

  if (D) {
    STAT(upolynomial, gcd_heuristic_success)++;
    if (debug_trace_ops.is_enabled("gcd")) {
      tracef("upolynomial_gcd_heuristic(%P, %P) = %P\n", A, B, D);
    }
  } else {
    if (debug_trace_ops.is_enabled("gcd")) {
      tracef("upolynomial_gcd_heuristic(%P, %P) failed", A, B, D);
    }
  }

  integer_ops.destruct(&A_cont);
  integer_ops.destruct(&B_cont);
  integer_ops.destruct(&d);
  integer_ops.destruct(&A_v);
  integer_ops.destruct(&B_v);
  integer_ops.destruct(&D_v);

  return D;
}
