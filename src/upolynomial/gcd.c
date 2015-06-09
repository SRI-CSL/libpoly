/*
 * univariate_polynomial_gcd.c
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#include "upolynomial/gcd.h"
#include "upolynomial/upolynomial_dense.h"
#include "upolynomial/upolynomial.h"
#include "upolynomial/output.h"

#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <assert.h>

STAT_DECLARE(int, upolynomial, gcd_euclid)
STAT_DECLARE(int, upolynomial, gcd_euclid_extended)
STAT_DECLARE(int, upolynomial, gcd_subresultant)
STAT_DECLARE(int, upolynomial, gcd_heuristic)
STAT_DECLARE(int, upolynomial, gcd_heuristic_success)

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
lp_upolynomial_t* upolynomial_gcd_euclid(const lp_upolynomial_t* A, const lp_upolynomial_t* B, lp_upolynomial_t** U, lp_upolynomial_t** V)
{
  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd_euclid("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(")\n");
  }

  assert(!lp_upolynomial_is_zero(B));

  int extended_gcd = (U != 0 && V != 0);

  if (extended_gcd) STAT(upolynomial, gcd_euclid_extended) ++;
  else STAT(upolynomial, gcd_euclid) ++;

  // Degrees of p and q
  size_t deg_A = lp_upolynomial_degree(A);
  assert(deg_A >= lp_upolynomial_degree(B));

  // The ring of computation
  assert(A->K == B->K);
  lp_int_ring_t* K = A->K;
  assert(K && K->is_prime);

  // The remainder in the calculation
  lp_upolynomial_t* D = 0;

  // Buffers to keep A, B, and remainders
  upolynomial_dense_t r_0;
  upolynomial_dense_t r_1;
  upolynomial_dense_construct_p(&r_0, deg_A + 1, A);
  upolynomial_dense_construct_p(&r_1, deg_A + 1, B);

  // Buffers to keep div and rem
  upolynomial_dense_t r_2;
  upolynomial_dense_t q;
  upolynomial_dense_construct(&r_2, deg_A + 1); // Some extra for reuse
  upolynomial_dense_construct(&q, deg_A + 1);

  // Buffers for the extended gcd computation
  upolynomial_dense_t s_0, s_1;
  upolynomial_dense_t t_0, t_1;
  if (extended_gcd) {
    upolynomial_dense_construct(&s_0, deg_A + 1);
    upolynomial_dense_construct(&s_1, deg_A + 1);
    upolynomial_dense_construct(&t_0, deg_A + 1);
    upolynomial_dense_construct(&t_1, deg_A + 1);
    // s_0, t_1 = 1
    integer_assign_int(lp_Z, s_0.coefficients, 1);
    integer_assign_int(lp_Z, t_1.coefficients, 1);
  }

  do {
    // One step of division
    upolynomial_dense_div_general(K, 1 /* exact */, &r_0, &r_1, &q, &r_2);

    if (trace_is_enabled("gcd")) {
      tracef("r_0 = ");
      upolynomial_dense_print(&r_0, trace_out);
      tracef("\nr_1 = ");
      upolynomial_dense_print(&r_1, trace_out);
      tracef("\nq = ");
      upolynomial_dense_print(&q, trace_out);
      tracef("\nr_2 = ");
      upolynomial_dense_print(&r_2, trace_out);
      if (extended_gcd) {
        tracef("\ns_0 = ");
        upolynomial_dense_print(&s_0, trace_out);
        tracef("\ns_1 = ");
        upolynomial_dense_print(&s_1, trace_out);
        tracef("\nt_0 = ");
        upolynomial_dense_print(&t_0, trace_out);
        tracef("\nt_1 = ");
        upolynomial_dense_print(&t_1, trace_out);
      }
      tracef("\n");
    }

    // Check we are done
    if (upolynomial_dense_is_zero(&r_2))  {
      // We're in a field, make it monic
      lp_integer_t lc;
      integer_construct_copy(K, &lc, r_1.coefficients + r_1.size - 1);
      if (integer_cmp_int(lp_Z, &lc, 1)) {
        upolynomial_dense_div_c(&r_1, K, &lc);
        if (extended_gcd) {
          upolynomial_dense_div_c(&s_1, K, &lc);
          upolynomial_dense_div_c(&t_1, K, &lc);
        }
      }
      integer_destruct(&lc);
      D = upolynomial_dense_to_upolynomial(&r_1, K);
      if (extended_gcd) {
        *U = upolynomial_dense_to_upolynomial(&s_1, K);
        *V = upolynomial_dense_to_upolynomial(&t_1, K);
      }
    } else {

      // Extended gcd computation
      if (extended_gcd) {
        // s2 = s0 - q*s1
        // (s0, s1) = (s1, s2)
        upolynomial_dense_sub_mult(&s_0, K, &q, &s_1);
        upolynomial_dense_swap(&s_0, &s_1);;
        // t2 = t0 - q*t2
        upolynomial_dense_sub_mult(&t_0, K, &q, &t_1);
        upolynomial_dense_swap(&t_0, &t_1);
      }

      // (r0, r1) = (r1, r2)
      upolynomial_dense_swap(&r_0, &r_1);
      upolynomial_dense_swap(&r_1, &r_2);


    }

  } while (D == 0);

  if (extended_gcd) {
    upolynomial_dense_destruct(&s_0);
    upolynomial_dense_destruct(&s_1);
    upolynomial_dense_destruct(&t_0);
    upolynomial_dense_destruct(&t_1);
  }

  upolynomial_dense_destruct(&q);
  upolynomial_dense_destruct(&r_2);
  upolynomial_dense_destruct(&r_0);
  upolynomial_dense_destruct(&r_1);

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd_euclid("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(") = "); lp_upolynomial_print(D, trace_out); tracef("\n");
  }

  return D;
}

lp_upolynomial_t* upolynomial_gcd_subresultant(const lp_upolynomial_t* A, const lp_upolynomial_t* B) {

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd_subresultant("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(")\n");
  }
  STAT(upolynomial, gcd_subresultant) ++;

  assert(!lp_upolynomial_is_zero(B));

  // The ring of compuation
  assert(A->K == B->K);
  assert(A->K == lp_Z);
  lp_int_ring_t* K = A->K;

  // Degrees of p and q
  size_t deg_A = lp_upolynomial_degree(A);
  assert(deg_A >= lp_upolynomial_degree(B));

  // The remainder in the calculation
  lp_upolynomial_t* D = 0;

  // Dense representations of the remainders to keep p and q
  upolynomial_dense_t r_0;
  upolynomial_dense_t r_1;
  upolynomial_dense_construct_p(&r_0, deg_A + 1, A);
  upolynomial_dense_construct_p(&r_1, deg_A + 1, B);

  // Contents of a, b
  lp_integer_t A_cont, B_cont;
  integer_construct_from_int(K, &A_cont, 0);
  integer_construct_from_int(K, &B_cont, 0);
  lp_upolynomial_content_Z(A, &A_cont);
  lp_upolynomial_content_Z(B, &B_cont);

  if (trace_is_enabled("gcd")) {
    tracef("cont(p) = "); integer_print(&A_cont, trace_out); tracef("\n");
    tracef("cont(q) = "); integer_print(&B_cont, trace_out); tracef("\n");
  }

  // d = gcd(content(p), content(q)))
  lp_integer_t d;
  integer_construct_from_int(K, &d, 1);

  integer_gcd_Z(&d, &A_cont, &B_cont);
  if (integer_cmp_int(lp_Z, &d, 1)) {
    // GCD != 1
    upolynomial_dense_div_c(&r_0, K, &A_cont);
    upolynomial_dense_div_c(&r_1, K, &B_cont);
  }

  // Buffers to keep div and rem
  upolynomial_dense_t r_2;
  upolynomial_dense_t q;
  upolynomial_dense_construct(&r_2, deg_A + 1); // Some extra for reuse
  upolynomial_dense_construct(&q, deg_A + 1);

  // Adjustment coefficients
  lp_integer_t g, h;
  integer_construct_from_int(lp_Z, &g, 1);
  integer_construct_from_int(lp_Z, &h, 1);

  // Temps for computation
  lp_integer_t tmp1, tmp2;
  integer_construct_from_int(K, &tmp1, 0);
  integer_construct_from_int(K, &tmp2, 0);

  do {
    // \delta = deg(p) - deg(q)
    int delta = r_0.size - r_1.size;

    // One step of division
    upolynomial_dense_div_general(K, 0, &r_0, &r_1, &q, &r_2);

    if (trace_is_enabled("gcd")) {
      tracef("r_0 = ");
      upolynomial_dense_print(&r_0, trace_out);
      tracef("\nr_q = ");
      upolynomial_dense_print(&r_1, trace_out);
      tracef("\nq = ");
      upolynomial_dense_print(&q, trace_out);
      tracef("\nr_w = ");
      upolynomial_dense_print(&r_2, trace_out);
      tracef("\n");
    }

    // Check if the remainder is of degree 0
    if (r_2.size == 1)  {
      if (integer_sgn(lp_Z, r_2.coefficients)) {
        // rem != 0, GCD(p, q) is 1 => total gcd is d
        integer_assign(K, r_2.coefficients, &d);
        D = upolynomial_dense_to_upolynomial(&r_2, K);
      } else {
        // rem == 0
        upolynomial_dense_mk_primitive_Z(&r_1, 1);
        upolynomial_dense_mult_c(&r_1, K, &d);
        D = upolynomial_dense_to_upolynomial(&r_1, K);
      }
    } else {
      // p = q
      upolynomial_dense_swap(&r_0, &r_1);
      // q = rem/(g*(h^delta)
      integer_pow(K, &tmp1, &h, delta);
      integer_mul(K, &tmp2, &tmp1, &g);
      upolynomial_dense_swap(&r_1, &r_2);
      upolynomial_dense_div_c(&r_1, K, &tmp2);
      // g = lc(p)
      integer_assign(K, &g, r_0.coefficients + r_0.size - 1);
      // h = h^(1-delta)*g^delta = g^delta/(h^(delta-1))
      integer_pow(K, &tmp1, &g, delta);
      integer_pow(K, &tmp2, &h, delta - 1);
      integer_div_exact(K, &h, &tmp1, &tmp2);
    }
  } while (D == 0);

  integer_destruct(&tmp1);
  integer_destruct(&tmp2);
  integer_destruct(&g);
  integer_destruct(&h);

  integer_destruct(&d);
  integer_destruct(&A_cont);
  integer_destruct(&B_cont);

  upolynomial_dense_destruct(&q);
  upolynomial_dense_destruct(&r_2);
  upolynomial_dense_destruct(&r_0);
  upolynomial_dense_destruct(&r_1);

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd_subresultant("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(") = "); lp_upolynomial_print(D, trace_out); tracef("\n");
  }

  return D;
}

static void evaluate_polynomial(const lp_upolynomial_t* A, const lp_integer_t* A_content, unsigned pow, lp_integer_t* out) {

  assert(pow > 0);

  integer_assign_int(lp_Z, out, 0);

  lp_integer_t add;
  lp_integer_t coeff;
  integer_construct_from_int(lp_Z, &add, 0);
  integer_construct_from_int(lp_Z, &coeff, 0);

  size_t k;
  for (k = 0; k < A->size; ++ k) {
    integer_div_exact(lp_Z, &coeff, &A->monomials[k].coefficient, A_content);
    if (A->monomials[k].degree == 0) {
      integer_assign(lp_Z, &add, &coeff);
    } else {
      integer_mul_pow2(lp_Z, &add, &coeff, A->monomials[k].degree * pow);
    }
    integer_add(lp_Z, out, out, &add);

    if (trace_is_enabled("gcd")) {
      tracef("out = "); integer_print(out, trace_out); tracef("\n");
    }
  }

  integer_destruct(&coeff);
  integer_destruct(&add);
}

/**
 * Reconstruct the polynomial from the value at 2^power, while multiplying with cont.
 */
static lp_upolynomial_t* reconstruct_polynomial(size_t max_size, lp_integer_t* p_value, unsigned pow, const lp_integer_t* cont) {

  assert(pow > 0);

  // computation here
  upolynomial_dense_t p_d;
  upolynomial_dense_construct(&p_d, max_size);

  lp_integer_t div;
  lp_integer_t rem;
  integer_construct_from_int(lp_Z, &div, 0);
  integer_construct_from_int(lp_Z, &rem, 0);

  lp_integer_t P;
  integer_construct_from_int(lp_Z, &P, 2);
  integer_pow(lp_Z, &P, &P, pow);

  // Basically compute digits modulo x_value
  int d = 0;
  while (integer_sgn(lp_Z, p_value)) {
    integer_div_rem_pow2_Z(&div, &rem, p_value, pow);

    if (integer_bits(&rem) + 1 >= pow) {
      // If biger than 2^n/2 then take subtract 2^n-1
      integer_sub(lp_Z, &rem, &rem, &P);
      integer_inc(lp_Z, &div);
    }

    integer_swap(&rem, p_d.coefficients + d);
    integer_swap(&div, p_value);

    d ++;
  }
  upolynomial_dense_touch(&p_d, d - 1);

  // We only care about primitive GCDs in the reconstruction
  upolynomial_dense_mk_primitive_Z(&p_d, 1);

  // Now, multiply with cont
  upolynomial_dense_mult_c(&p_d, lp_Z, cont);

  // Get the sparse representation
  lp_upolynomial_t* p = upolynomial_dense_to_upolynomial(&p_d, lp_Z);

  // Free temporaries
  upolynomial_dense_destruct(&p_d);
  integer_destruct(&div);
  integer_destruct(&rem);
  integer_destruct(&P);

  return p;
}

int bound_valuation(const lp_upolynomial_t* A, const lp_upolynomial_t* B, const lp_integer_t* A_cont, const lp_integer_t* B_cont) {
  // Get the highest power 2^n such that
  //  2^n >= 2*min(|pp(A)|, |pp(B)|) + 2

  int A_max = 0;
  int B_max = 0;
  size_t k;

  lp_integer_t tmp;
  integer_construct_from_int(lp_Z, &tmp, 0);

  for (k = 0; k < A->size; ++ k) {
    integer_div_Z(&tmp, &A->monomials[k].coefficient, A_cont);
    int bits = integer_bits(&tmp);
    if (bits > A_max) {
      A_max = bits;
    }
  }

  for (k = 0; k < B->size; ++ k) {
    integer_div_Z(&tmp, &B->monomials[k].coefficient, B_cont);
    int bits = integer_bits(&tmp);
    if (bits > B_max) {
      B_max = bits;
    }
  }

  integer_destruct(&tmp);

  int bits_min = A_max > B_max ? B_max : A_max;

  return bits_min + 2;
}

lp_upolynomial_t* upolynomial_gcd_heuristic(const lp_upolynomial_t* A, const lp_upolynomial_t* B, int attempts) {

  // Let's keep the smaller one in B
  if (lp_upolynomial_degree(A) < lp_upolynomial_degree(B)) {
    return upolynomial_gcd_heuristic(B, A, attempts);
  }

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd_heuristic("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(")\n");
  }
  STAT(upolynomial, gcd_heuristic) ++;

  lp_upolynomial_t* D = 0;

  // The ring of computation
  assert(A->K == B->K);
  assert(A->K == lp_Z);

  // content(A), content(B)
  lp_integer_t A_cont, B_cont;
  integer_construct_from_int(lp_Z, &A_cont, 0);
  integer_construct_from_int(lp_Z, &B_cont, 0);
  lp_upolynomial_content_Z(A, &A_cont);
  lp_upolynomial_content_Z(B, &B_cont);

  if (trace_is_enabled("gcd")) {
    tracef("cont(p) = "); integer_print(&A_cont, trace_out); tracef("\n");
    tracef("cont(q) = "); integer_print(&B_cont, trace_out); tracef("\n");
  }

  // d = gcd(content(A), content(B)))
  lp_integer_t d;
  integer_construct_from_int(lp_Z, &d, 1);
  integer_gcd_Z(&d, &A_cont, &B_cont);

  if (trace_is_enabled("gcd")) {
    tracef("d = "); integer_print(&d, trace_out); tracef("\n");
  }

  // The number we use for valuation 2^n
  int n = bound_valuation(A, B, &A_cont, &B_cont);

  lp_integer_t A_v, B_v, D_v;
  integer_construct_from_int(lp_Z, &A_v, 0);
  integer_construct_from_int(lp_Z, &B_v, 0);
  integer_construct_from_int(lp_Z, &D_v, 0);

  while (D == 0 && (attempts --)) {

    // Evaluate A and B
    evaluate_polynomial(A, &A_cont, n, &A_v);
    evaluate_polynomial(B, &B_cont, n, &B_v);

    if (trace_is_enabled("gcd")) {
        tracef("value of A/cont = "); integer_print(&A_v, trace_out); tracef("\n");
        tracef("value of B/cont = "); integer_print(&B_v, trace_out); tracef("\n");
    }

    // Get the gcd of the values and reconstruct the possible gcd
    integer_gcd_Z(&D_v, &A_v, &B_v);
    // This also changes the value D_v (B is the smaller one, but we need the size of A for reconstruction)
    D = reconstruct_polynomial(lp_upolynomial_degree(A) + 1, &D_v, n, &d);

    // Check if it divides both (B is smaller degree)
    if (!lp_upolynomial_divides(D, B) || !lp_upolynomial_divides(D, A)) {
      lp_upolynomial_delete(D);
      D = 0;
    }

    // Try next power
    n ++;
  }

  if (D) {
    STAT(upolynomial, gcd_heuristic_success)++;
    if (trace_is_enabled("gcd")) {
      tracef("upolynomial_gcd_heuristic("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(") = "); lp_upolynomial_print(D, trace_out); tracef("\n");
    }
  } else {
    if (trace_is_enabled("gcd")) {
      tracef("upolynomial_gcd_heuristic("); lp_upolynomial_print(A, trace_out); tracef(", "); lp_upolynomial_print(B, trace_out); tracef(") failed");
    }
  }

  integer_destruct(&A_cont);
  integer_destruct(&B_cont);
  integer_destruct(&d);
  integer_destruct(&A_v);
  integer_destruct(&B_v);
  integer_destruct(&D_v);

  return D;
}
