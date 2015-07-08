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
void upolynomial_root_bound_cauchy(const lp_upolynomial_t* f, lp_integer_t* B) {

  assert(f->K == lp_Z);

  lp_integer_t M_p; // Positive
  lp_integer_t M_n; // Negative
  integer_construct_from_int(lp_Z, &M_p, 0);
  integer_construct_from_int(lp_Z, &M_n, 0);

  int k;
  int d = f->size - 1;
  for (k = 0; k < d; ++ k) {
    if (f->monomials[k].degree > 0) {
      if (integer_cmp(lp_Z, &f->monomials[k].coefficient, &M_p) > 0) {
        integer_assign(lp_Z, &M_p, &f->monomials[k].coefficient);
        integer_neg(lp_Z, &M_n, &f->monomials[k].coefficient);
      } else if (integer_cmp(lp_Z, &f->monomials[k].coefficient, &M_n) < 0) {
        integer_assign(lp_Z, &M_n, &f->monomials[k].coefficient);
        integer_neg(lp_Z, &M_p, &f->monomials[k].coefficient);
      }
    }
  }

  integer_construct_from_int(lp_Z, B, 0);
  if (integer_sgn(lp_Z, &f->monomials[d].coefficient) > 0) {
    integer_div_Z(B, &M_p, &f->monomials[d].coefficient);
  } else {
    lp_integer_t neg;
    integer_construct_from_int(lp_Z, &neg, 0);
    integer_neg(lp_Z, &neg, &f->monomials[d].coefficient);
    integer_div_Z(B, &M_p, &f->monomials[d].coefficient);
    integer_destruct(&neg);
  }

  integer_inc(lp_Z, B);

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
void upolynomial_factor_bound_landau_mignotte(const lp_upolynomial_t* f, size_t n, lp_integer_t* B) {

  assert(f->K == lp_Z);
  assert(lp_upolynomial_degree(f) >= n);

  lp_integer_t tmp;
  lp_integer_t norm;
  lp_integer_t ak_sq;
  integer_construct_from_int(lp_Z, &tmp, 0);
  integer_construct_from_int(lp_Z, &norm, 0);
  integer_construct_from_int(lp_Z, &ak_sq, 0);

  size_t k;

  // Sum of squares
  for(k = 0; k < f->size; ++ k) {
    const lp_integer_t* ak = &f->monomials[k].coefficient;
    integer_mul(lp_Z, &ak_sq, ak, ak);
    integer_add(lp_Z, &tmp, &norm, &ak_sq);
    integer_swap(&tmp, &norm);
  }

  // Square root
  integer_sqrt_Z(&tmp, &norm);
  integer_swap(&tmp, &norm);
  integer_inc(lp_Z, &norm);

  // Power of two
  lp_integer_t two;
  lp_integer_t power;
  integer_construct_from_int(lp_Z, &two, 2);
  integer_construct_from_int(lp_Z, &power, 0);
  integer_pow(lp_Z, &power, &two, n);

  // Result
  integer_mul(lp_Z, B, &power, &norm);

  integer_destruct(&tmp);
  integer_destruct(&norm);
  integer_destruct(&ak_sq);
  integer_destruct(&two);
  integer_destruct(&power);
}
