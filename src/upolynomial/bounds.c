/*
 * univariate_polynomial_bounds.c
 *
 *  Created on: Dec 6, 2013
 *      Author: dejan
 */

#include "upolynomial/bounds.h"

#include <assert.h>

/**
 * Let f = a_0 + ... + a_n x^n and let z be a root.
 *
 * If |z| <= 1 the bound holds.
 *
 * If |z| > 1 then
 *
 *   a_n z^n = -a_0 - ... - a_1*z^{n-1}
 *
 *   |a_n| |z|^n <= max(a_0, ..., a_{n-1}) (1 + ... + |z|^{n-1})
 *                < max * |z|^n / (|z| - 1)
 *
 *   |z| < (max / |a_n|) + 1
 */
void upolynomial_root_bound_cauchy(const upolynomial_t* f, integer_t* B) {

  assert(f->K == Z);

  integer_t M_p; // Positive
  integer_t M_n; // Negative
  integer_construct_from_int(Z, &M_p, 0);
  integer_construct_from_int(Z, &M_n, 0);

  int k;
  int d = f->size - 1;
  for (k = 0; k < d; ++ k) {
    if (f->monomials[k].degree > 0) {
      if (integer_cmp(Z, &f->monomials[k].coefficient, &M_p) > 0) {
        integer_assign(Z, &M_p, &f->monomials[k].coefficient);
        integer_neg(Z, &M_n, &f->monomials[k].coefficient);
      } else if (integer_cmp(Z, &f->monomials[k].coefficient, &M_n) < 0) {
        integer_assign(Z, &M_n, &f->monomials[k].coefficient);
        integer_neg(Z, &M_p, &f->monomials[k].coefficient);
      }
    }
  }

  integer_construct_from_int(Z, B, 0);
  if (integer_sgn(Z, &f->monomials[d].coefficient) > 0) {
    integer_div_Z(B, &M_p, &f->monomials[d].coefficient);
  } else {
    integer_t neg;
    integer_construct_from_int(Z, &neg, 0);
    integer_neg(Z, &neg, &f->monomials[d].coefficient);
    integer_div_Z(B, &M_p, &f->monomials[d].coefficient);
    integer_destruct(&neg);
  }

  integer_inc(Z, B);

  integer_destruct(&M_p);
  integer_destruct(&M_n);
}

/**
 * Let
 *
 *  f = a_n*x^m + ... + a_0
 *  g = b_n*x^n + ... + b_0
 *
 * and g divides f.
 *
 * Landau's inequality: For any f, let z_1, ..., z_m be the (complex) roots of
 * f, then
 *
 *   M(f) = |a_n|*(max(1, |z_1|)*...*max(1,|z_m|)) <= norm(f).
 *
 * Landau proof: First, for any polynomial h, and any complex z it holds that
 *
 * [1] norm((x + z)*h) = norm(conj(z)x + 1)h).
 *
 * This can be checked by simple expansion. Now, let z_1, ..., z_k be the roots
 * with |z_i| > 1 and set
 *
 *   f' = a_n * (conj(z_1)*x - 1) * ... * (conj(z_k)x - 1)
 *            * (x - z_{k+1)      * ... * (x - z_m)
 * Then applying [1] we get that
 *
 *   norm(f)^2 = norm(f')^2 >= |lc(f')|^2 = M(f)^2.
 *
 * Bound proof: First notice that
 *
 * [2] |b_0| + ... + |b_n| <= 2^n*M(g).
 *
 * To see this, look at g as b_n*(x - z_1)*...*(x - z_n) and therefore
 *
 *  |b_j| = |b_n|*|\sum x^j*z_i1*...*z_i(n-j)| <= binomial(n, j)*M(g)
 *
 * Summing up above we get [2].
 *
 * Since g divides f we also know that
 *
 *  M(g) <= |b_n/a_m|*M(f) <= |b_n/a_m|*norm(f) (by Landau)
 *
 * Now we have enouth to get
 *
 *  |b_j| <= 2^n*norm(f).
 */
void upolynomial_factor_bound_landau_mignotte(const upolynomial_t* f, size_t n, integer_t* B) {

  assert(f->K == Z);
  assert(upolynomial_ops.degree(f) >= n);

  integer_t tmp;
  integer_t norm;
  integer_t ak_sq;
  integer_construct_from_int(Z, &tmp, 0);
  integer_construct_from_int(Z, &norm, 0);
  integer_construct_from_int(Z, &ak_sq, 0);

  int k;

  // Sum of squares
  for(k = 0; k < f->size; ++ k) {
    const integer_t* ak = &f->monomials[k].coefficient;
    integer_mul(Z, &ak_sq, ak, ak);
    integer_add(Z, &tmp, &norm, &ak_sq);
    integer_swap(Z, &tmp, &norm);
  }

  // Square root
  integer_sqrt_Z(&tmp, &norm);
  integer_swap(Z, &tmp, &norm);
  integer_inc(Z, &norm);

  // Power of two
  integer_t two;
  integer_t power;
  integer_construct_from_int(Z, &two, 2);
  integer_construct_from_int(Z, &power, 0);
  integer_pow(Z, &power, &two, n);

  // Result
  integer_mul(Z, B, &power, &norm);

  integer_destruct(&tmp);
  integer_destruct(&norm);
  integer_destruct(&ak_sq);
  integer_destruct(&two);
  integer_destruct(&power);
}
