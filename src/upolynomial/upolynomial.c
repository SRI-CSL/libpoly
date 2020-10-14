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

#include "upolynomial/umonomial.h"
#include "upolynomial/upolynomial_dense.h"

#include "upolynomial/output.h"
#include "upolynomial/gcd.h"
#include "upolynomial/factorization.h"
#include "upolynomial/root_finding.h"

#include "utils/debug_trace.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <upolynomial.h>
#include "upolynomial/upolynomial.h"

size_t lp_upolynomial_degree(const lp_upolynomial_t* p) {
  assert(p);
  assert(p->size > 0);
  return p->monomials[p->size-1].degree;
}

// Construct the polynomial, but don't construct monomials
lp_upolynomial_t* lp_upolynomial_construct_empty(const lp_int_ring_t* K, size_t size) {
  size_t malloc_size = sizeof(lp_upolynomial_t) + size*sizeof(ulp_monomial_t);
  lp_upolynomial_t* new_p = (lp_upolynomial_t*) malloc(malloc_size);
  new_p->K = (lp_int_ring_t*)K;
  new_p->size = size;
  lp_int_ring_attach((lp_int_ring_t*)K);
  return new_p;
}

lp_upolynomial_t* lp_upolynomial_construct(const lp_int_ring_t* K, size_t degree, const lp_integer_t* coefficients) {

  // Compute the needed size
  unsigned i;
  unsigned size = 0;

  for (i = 0; i <= degree; ++ i) {
    if (integer_sgn(K, coefficients + i)) {
      size ++;
    }
  }

  if (size == 0) {
    size = 1;
  }

  // Allocate the memory (at least one)
  lp_upolynomial_t* new_p = lp_upolynomial_construct_empty(K, size > 0 ? size : 1);

  // Copy the coefficients
  size = 0;
  lp_integer_t tmp;
  integer_construct_from_int(lp_Z, &tmp, 0);
  for (i = 0; i <= degree; ++ i) {
    integer_assign(K, &tmp, coefficients + i);
    if (integer_sgn(lp_Z, &tmp)) {
      umonomial_construct(K, &new_p->monomials[size], i, &tmp);
      size ++;
    }
  }
  integer_destruct(&tmp);

  // If no one constructed, make 0
  if (size == 0) {
    umonomial_construct_from_int(K, &new_p->monomials[0], 0, 0);
    size = 1;
  }

  // Return the new polynomial
  return new_p;
}

void lp_upolynomial_delete(lp_upolynomial_t* p) {
  assert(p);
  size_t i = 0;
  for (i = 0; i < p->size; ++ i) {
    integer_destruct(&p->monomials[i].coefficient);
  }
  lp_int_ring_detach((lp_int_ring_t*)p->K);
  free(p);
}

lp_upolynomial_t* lp_upolynomial_construct_power(const lp_int_ring_t* K, size_t degree, long c) {
  lp_upolynomial_t* result = lp_upolynomial_construct_empty(K, 1);
  integer_construct_from_int(K, &result->monomials[0].coefficient, c);
  result->monomials[0].degree = degree;
  return result;
}

lp_upolynomial_t* lp_upolynomial_construct_from_int(const lp_int_ring_t* K, size_t degree, const int* coefficients) {

  unsigned i;
  lp_integer_t real_coefficients[degree+1];

  for (i = 0; i <= degree; ++ i) {
    integer_construct_from_int(K, &real_coefficients[i], coefficients[i]);
  }

  lp_upolynomial_t* result = lp_upolynomial_construct(K, degree, real_coefficients);

  for (i = 0; i <= degree; ++ i) {
    integer_destruct(&real_coefficients[i]);
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_construct_from_long(const lp_int_ring_t* K, size_t degree, const long* coefficients) {

  unsigned i;
  lp_integer_t real_coefficients[degree+1];

  for (i = 0; i <= degree; ++ i) {
    integer_construct_from_int(K, &real_coefficients[i], coefficients[i]);
  }

  lp_upolynomial_t* result = lp_upolynomial_construct(K, degree, real_coefficients);

  for (i = 0; i <= degree; ++ i) {
    integer_destruct(&real_coefficients[i]);
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_construct_copy(const lp_upolynomial_t* p) {
  assert(p);
  // Allocate the memory
  lp_upolynomial_t* new_p = lp_upolynomial_construct_empty(p->K, p->size);
  // Copy the coefficients
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    umonomial_construct_copy(lp_Z, &new_p->monomials[i], &p->monomials[i]);
  }
  // Return the new polynomial
  return new_p;
}

lp_upolynomial_t* lp_upolynomial_construct_copy_K(const lp_int_ring_t* K, const lp_upolynomial_t* p) {

  assert(p);
  assert(K != p->K);

  // If in Z, same, or bigger, it's just a regular copy (while swapping the rings).
  if (K == lp_Z || (p->K != lp_Z && integer_cmp(lp_Z, &K->M, &p->K->M) >= 0)) {
    lp_upolynomial_t* copy = lp_upolynomial_construct_copy(p);
    lp_int_ring_detach(copy->K);
    copy->K = (lp_int_ring_t*)K;
    lp_int_ring_attach(copy->K);
    return copy;
  }

  // Compute the needed size
  unsigned i;
  size_t size = 0;
  size_t degree = 0;

  for (i = 0; i < p->size; ++ i) {
    if (integer_sgn(K, &p->monomials[i].coefficient)) {
      size ++;
      degree = p->monomials[i].degree;
    }
  }

  // At least one element
  if (degree == 0) {
    size = 1;
  }

  // Allocate the memory
  lp_upolynomial_t* new_p = lp_upolynomial_construct_empty(K, size);

  // Copy the coefficients
  if (degree == 0) {
    if (integer_sgn(K, &p->monomials[0].coefficient)) {
      umonomial_construct_copy(K, &new_p->monomials[0], &p->monomials[0]);
    } else {
      umonomial_construct_from_int(K, &new_p->monomials[0], 0, 0);
    }
  } else {
    size = 0;
    lp_integer_t tmp;
    integer_construct_from_int(lp_Z, &tmp, 0);
    for (i = 0; i < p->size; ++ i) {
      integer_assign(K, &tmp, &p->monomials[i].coefficient);
      if (integer_sgn(lp_Z, &tmp)) {
        umonomial_construct(lp_Z, &new_p->monomials[size], p->monomials[i].degree, &tmp);
        size ++;
      }
    }
    integer_destruct(&tmp);
  }

  // Return the new polynomial
  return new_p;
}

void lp_upolynomial_unpack(const lp_upolynomial_t* p, lp_integer_t* out) {
  assert(p);
  unsigned i; // i big step, j small step
  for (i = 0; i < p->size; ++ i) {
    integer_assign(lp_Z, &out[p->monomials[i].degree], &p->monomials[i].coefficient);
  }
}

const lp_int_ring_t* lp_upolynomial_ring(const lp_upolynomial_t* p) {
  assert(p);
  return p->K;
}

void lp_upolynomial_set_ring(lp_upolynomial_t* p, const lp_int_ring_t* K) {
  assert(p);
  lp_int_ring_detach(p->K);
  p->K = (lp_int_ring_t*)K;
  lp_int_ring_attach(p->K);
}

const lp_integer_t* lp_upolynomial_lead_coeff(const lp_upolynomial_t* p) {
  assert(p);
  assert(p->size > 0);
  return &p->monomials[p->size-1].coefficient;
}

const lp_integer_t* lp_upolynomial_const_term(const lp_upolynomial_t* p) {
  assert(p);
  if (p->monomials[0].degree == 0) {
    return &p->monomials[0].coefficient;
  } else {
    return 0;
  }
}

int lp_upolynomial_cmp(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {
  assert(p);
  assert(q);
  assert(p->K == q->K);

  int p_deg, q_deg, p_i, q_i, cmp;

  p_i = p->size;
  q_i = q->size;

  do {

    // Next coefficient to compare
    p_i --;
    q_i --;

    // Compare degrees
    p_deg = p->monomials[p_i].degree;
    q_deg = q->monomials[q_i].degree;
    if (p_deg > q_deg) return 1;
    if (q_deg > p_deg) return -1;

    // Degrees equal, compare coefficients
    cmp = integer_cmp(lp_Z, &p->monomials[p_i].coefficient, &q->monomials[q_i].coefficient);
    if (cmp) {
      return cmp;
    }

    // Both equal, go to the next one if there is one
  } while (p_i > 0 || q_i > 0);

  // One of the polynomials ran out of coefficients
  if (p_i == q_i) return 0;
  else if (p_i > q_i) return 1;
  else return -1;
}

int lp_upolynomial_is_zero(const lp_upolynomial_t* p) {
  if (p->size > 1) return 0;
  if (p->monomials[0].degree > 0) return 0;
  return integer_sgn(lp_Z, &p->monomials[0].coefficient) == 0;
}

int lp_upolynomial_is_one(const lp_upolynomial_t* p) {
  if (p->size > 1) return 0;
  if (p->monomials[0].degree > 0) return 0;
  return integer_cmp_int(lp_Z, &p->monomials[0].coefficient, 1) == 0;
}

int lp_upolynomial_is_monic(const lp_upolynomial_t* p) {
  return integer_cmp_int(lp_Z, lp_upolynomial_lead_coeff(p), 1) == 0;
}

lp_upolynomial_t* lp_upolynomial_add(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_add("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  const lp_int_ring_t* K = p->K;

  // Get the end degree
  size_t degree = lp_upolynomial_degree(p);
  if (lp_upolynomial_degree(q) > degree) {
    degree = lp_upolynomial_degree(q);
  }

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_construct(&tmp, degree + 1);

  // Add the buffers
  upolynomial_dense_add_mult_p_int(&tmp, p, 1);
  upolynomial_dense_add_mult_p_int(&tmp, q, 1);

  // Construct the result
  lp_upolynomial_t* result = upolynomial_dense_to_upolynomial(&tmp, K);
  upolynomial_dense_destruct(&tmp);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_add("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_sub(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_sub("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  const lp_int_ring_t* K = p->K;

  // Get the end degree
  size_t degree = lp_upolynomial_degree(p);
  if (lp_upolynomial_degree(q) > degree) {
    degree = lp_upolynomial_degree(q);
  }

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_construct(&tmp, degree + 1);

  // Add the buffers
  upolynomial_dense_add_mult_p_int(&tmp, p, 1);
  upolynomial_dense_add_mult_p_int(&tmp, q, -1);

  // Construct the result
  lp_upolynomial_t* result = upolynomial_dense_to_upolynomial(&tmp, K);
  upolynomial_dense_destruct(&tmp);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_sub("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_multiply_simple(const ulp_monomial_t* m, const lp_upolynomial_t* q) {

  assert(m);
  assert(q);

  lp_upolynomial_t* result = lp_upolynomial_construct_copy(q);

  size_t i;
  for (i = 0; i < result->size; ++ i) {
    integer_mul(q->K, &result->monomials[i].coefficient, &m->coefficient, &q->monomials[i].coefficient);
    result->monomials[i].degree += m->degree;
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_mul(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_multiply("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  // Take p to be the smaller size
  if (p->size > q->size) {
    return lp_upolynomial_mul(q, p);
  }

  // If multiplying with zero
  if (lp_upolynomial_is_zero(p) || lp_upolynomial_is_zero(q)) {
    return lp_upolynomial_construct_power(p->K, 0, 0);
  }

  // Special case if possible
  if (p->K == lp_Z && p->size == 1) {
    lp_upolynomial_t* result = lp_upolynomial_multiply_simple(p->monomials, q);
    if (trace_is_enabled("arithmetic")) {
      tracef("upolynomial_multiply("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
    }
    return result;
  }

  // Max degree of the multiplications
  size_t degree = lp_upolynomial_degree(p)+lp_upolynomial_degree(q);

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_construct(&tmp, degree + 1);

  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    upolynomial_dense_add_mult_p_mon(&tmp, q, &p->monomials[i]);
  }

  // Construct the result
  lp_upolynomial_t* result = upolynomial_dense_to_upolynomial(&tmp, p->K);
  upolynomial_dense_destruct(&tmp);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_multiply("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_mul_c(const lp_upolynomial_t* p, const lp_integer_t* c) {
  assert(p);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_multiply_c(");
    lp_upolynomial_print(p, trace_out); tracef(", ");
    integer_print(c, trace_out); tracef(")\n");
  }

  ulp_monomial_t m;
  umonomial_construct(p->K, &m, 0, c);
  lp_upolynomial_t* result = lp_upolynomial_multiply_simple(&m, p);
  umonomial_destruct(&m);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_multiply_c(");
    lp_upolynomial_print(p, trace_out); tracef(", ");
    integer_print(c, trace_out); tracef(") = ");
    lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_pow(const lp_upolynomial_t* p, long pow) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_pow("); lp_upolynomial_print(p, trace_out); tracef(", %ld)\n", pow);
  }

  assert(p);
  assert(pow >= 0);

  lp_upolynomial_t* result  = 0;

  // Anything of size 1
  if (p->size == 1) {
    result = lp_upolynomial_construct_empty(p->K, 1);
    integer_construct_from_int(lp_Z, &result->monomials[0].coefficient, 0);
    integer_pow(p->K, &result->monomials[0].coefficient, &p->monomials[0].coefficient, pow);
    result->monomials[0].degree = p->monomials[0].degree * pow;
  } else {
    result = lp_upolynomial_construct_power(p->K, 0, 1);
    lp_upolynomial_t* tmp = lp_upolynomial_construct_copy(p);
    lp_upolynomial_t* prev;
    while (pow) {
      if (pow & 1) {
        prev = result;
        result = lp_upolynomial_mul(result, tmp);
        lp_upolynomial_delete(prev);
      }
      prev = tmp;
      tmp = lp_upolynomial_mul(tmp, tmp);
      pow >>= 1;
      lp_upolynomial_delete(prev);
    }
    lp_upolynomial_delete(tmp);
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_pow(");
    lp_upolynomial_print(p, trace_out);
    tracef(", %ld) = ", pow);
    lp_upolynomial_print(result, trace_out);
    tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_derivative(const lp_upolynomial_t* p) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_derivative("); lp_upolynomial_print(p, trace_out); tracef(")\n");
  }

  // Max degree of the derivative
  size_t degree = lp_upolynomial_degree(p);
  if (degree > 0) degree --;

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_construct(&tmp, degree + 1);
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    size_t degree_i = p->monomials[i].degree;
    if (degree_i > 0) {
      integer_mul_int(p->K, &tmp.coefficients[degree_i-1], &p->monomials[i].coefficient, degree_i);
    }
  }
  tmp.size = degree + 1;

  // Construct the result
  lp_upolynomial_t* result = upolynomial_dense_to_upolynomial(&tmp, p->K);
  upolynomial_dense_destruct(&tmp);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_derivative("); lp_upolynomial_print(p, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

// Performs general division loop, result is returned in the buffers
// Buffers should be passed in unconstructed
void lp_upolynomial_div_general(const lp_upolynomial_t* p, const lp_upolynomial_t* q, upolynomial_dense_t* div, upolynomial_dense_t* rem, int exact) {

  if (trace_is_enabled("division")) {
    tracef("upolynomial_div_general("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(lp_upolynomial_degree(q) <= lp_upolynomial_degree(p));

  const lp_int_ring_t* K = p->K;

  int p_deg = lp_upolynomial_degree(p);
  int q_deg = lp_upolynomial_degree(q);

  upolynomial_dense_construct_p(rem, p_deg + 1, p);
  upolynomial_dense_construct(div, p_deg - q_deg + 1);

  // monomial we use to multiply with
  ulp_monomial_t m;
  umonomial_construct_from_int(lp_Z, &m, 0, 0);

  // adjustment for the powers of div
  lp_integer_t adjust;
  integer_construct_from_int(lp_Z, &adjust, 0);

  int k;
  for (k = p_deg; k >= q_deg; -- k) {
    // next coefficient to eliminate
    while (exact && k >= q_deg && integer_sgn(lp_Z, rem->coefficients + k) == 0) {
      k --;
    }
    // if found, eliminate it
    if (k >= q_deg) {

      if (trace_is_enabled("division")) {
        tracef("dividing with "); lp_upolynomial_print(q, trace_out); tracef(" at degree %d\n", k);
        tracef("rem = "); upolynomial_dense_print(rem, trace_out);
        tracef("div = "); tracef("\n");
        upolynomial_dense_print(div, trace_out); tracef("\n");
      }

      assert(!exact || integer_divides(K, lp_upolynomial_lead_coeff(q), rem->coefficients + k));

      // Eliminate the coefficient using q*m
      m.degree = k - q_deg;

      // Get the quotient coefficient
      if (exact) {
        // get the quotient
        integer_div_exact(K, &m.coefficient, rem->coefficients + k, lp_upolynomial_lead_coeff(q));
      } else {
        // rem: a*x^k, q: b*x^d, so we multiply rem with b and subtract a*q
        integer_assign(lp_Z, &m.coefficient, rem->coefficients + k);
        upolynomial_dense_mult_c(rem, K, lp_upolynomial_lead_coeff(q));
      }

      // Do the subtraction
      if (integer_sgn(K, &m.coefficient)) {
        upolynomial_dense_sub_mult_p_mon(rem, q, &m);
      }

      // Put the monomial into the division
      if (exact || !integer_sgn(lp_Z, &m.coefficient)) {
        integer_swap(&div->coefficients[m.degree], &m.coefficient);
      } else {
        // Adjust the monomial with to be lc(q)^(k-deg_q)
        integer_pow(K, &adjust, lp_upolynomial_lead_coeff(q), m.degree);
        integer_mul(K, &div->coefficients[m.degree], &m.coefficient, &adjust);
      }
      upolynomial_dense_touch(div, m.degree);
    }
  }

  // destruct temps
  integer_destruct(&adjust);
  umonomial_destruct(&m);
}

lp_upolynomial_t* lp_upolynomial_div_degrees(const lp_upolynomial_t* p, size_t a) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_degrees("); lp_upolynomial_print(p, trace_out); tracef(", %zd)\n", a);
  }

  assert(a > 1);

  lp_upolynomial_t* result = lp_upolynomial_construct_copy(p);
  size_t i;
  for (i = 0; i < result->size; ++ i) {
    assert(result->monomials[i].degree % a == 0);
    result->monomials[i].degree /= a;
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_degrees("); lp_upolynomial_print(p, trace_out); tracef(", %zd) = ", a); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_div_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!lp_upolynomial_is_zero(q));

  lp_upolynomial_t* result = 0;

  if (lp_upolynomial_degree(p) >= lp_upolynomial_degree(q)) {
    const lp_int_ring_t* K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    lp_upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    result = upolynomial_dense_to_upolynomial(&div_buffer, K);
    upolynomial_dense_destruct(&div_buffer);
    upolynomial_dense_destruct(&rem_buffer);
  } else {
    // 0 polynomial
    result = lp_upolynomial_construct_power(p->K, 0, 0);
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_div_exact_c(const lp_upolynomial_t* p, const lp_integer_t* c) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact(");
    lp_upolynomial_print(p, trace_out); tracef(", ");
    integer_print(c, trace_out); tracef(")\n");
  }

  const lp_int_ring_t* K = p->K;

  assert(p);
  assert(integer_cmp_int(K, c, 0));

  lp_upolynomial_t* result = lp_upolynomial_construct_empty(K, p->size);

  size_t i;
  for (i = 0; i < p->size; ++ i) {
    result->monomials[i].degree = p->monomials[i].degree;
    integer_construct_from_int(K, &result->monomials[i].coefficient, 0);
    assert(integer_divides(K, c, &p->monomials[i].coefficient));
    integer_div_exact(K, &result->monomials[i].coefficient, &p->monomials[i].coefficient, c);
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact(");
    lp_upolynomial_print(p, trace_out); tracef(", ");
    integer_print(c, trace_out); tracef(") = ");
    lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

lp_upolynomial_t* lp_upolynomial_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!lp_upolynomial_is_zero(q));

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_rem_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  lp_upolynomial_t* result = 0;

  if (lp_upolynomial_degree(p) >= lp_upolynomial_degree(q)) {
    const lp_int_ring_t* K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    lp_upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    result = upolynomial_dense_to_upolynomial(&rem_buffer, K);
    upolynomial_dense_destruct(&rem_buffer);
    upolynomial_dense_destruct(&div_buffer);
  } else {
    // rem = p
    result = lp_upolynomial_construct_copy(p);
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_rem_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(result, trace_out); tracef("\n");
  }

  return result;
}

void lp_upolynomial_div_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q,
    lp_upolynomial_t** div, lp_upolynomial_t** rem) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_rem_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!lp_upolynomial_is_zero(q));
  assert(*div == 0 && *rem == 0);

  if (lp_upolynomial_degree(p) >= lp_upolynomial_degree(q)) {
    const lp_int_ring_t* K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    lp_upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    *div= upolynomial_dense_to_upolynomial(&div_buffer, K);
    *rem = upolynomial_dense_to_upolynomial(&rem_buffer, K);
    upolynomial_dense_destruct(&div_buffer);
    upolynomial_dense_destruct(&rem_buffer);
  } else {
    // 0 polynomial
    *div = lp_upolynomial_construct_power(p->K, 0, 0);
    // rem = p
    *rem = lp_upolynomial_construct_copy(p);
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = ("); lp_upolynomial_print(*div, trace_out); tracef(", "); lp_upolynomial_print(*rem, trace_out); tracef(")\n");
  }
}


void lp_upolynomial_div_pseudo(lp_upolynomial_t** div, lp_upolynomial_t** rem, const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_pseudo("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(!*div);
  assert(!*rem);
  assert(p->K == q->K);
  assert(!lp_upolynomial_is_zero(q));

  assert(lp_upolynomial_degree(p) >= lp_upolynomial_degree(q));

  const lp_int_ring_t* K = p->K;

  upolynomial_dense_t rem_buffer;
  upolynomial_dense_t div_buffer;

  lp_upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** pseudo */ 0);

  *div = upolynomial_dense_to_upolynomial(&div_buffer, K);
  *rem = upolynomial_dense_to_upolynomial(&rem_buffer, K);

  upolynomial_dense_destruct(&div_buffer);
  upolynomial_dense_destruct(&rem_buffer);

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_pseudo("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = ("); lp_upolynomial_print(*div, trace_out); tracef(", "); lp_upolynomial_print(*rem, trace_out); tracef(")\n");
  }
}

int lp_upolynomial_divides(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p->K == q->K);

  // Special case
  if (lp_upolynomial_degree(p) > lp_upolynomial_degree(q)) {
    return 0;
  }

  const lp_int_ring_t* K = p->K;

  int result = 0;

  // Check that the least power of x divides the the q
  if (p->monomials[0].degree > q->monomials[0].degree) {
    result = 0;
  } else {
    // Check that the least coefficient divides
    if (!integer_divides(K, &p->monomials[0].coefficient, &q->monomials[0].coefficient)) {
      result = 0;
    } else {
      // Let's just divide
      if (K && K->is_prime) {
        lp_upolynomial_t* rem = lp_upolynomial_rem_exact(q, p);
        result = lp_upolynomial_is_zero(rem);
        lp_upolynomial_delete(rem);
      } else {
        lp_upolynomial_t* div = 0;
        lp_upolynomial_t* rem = 0;
        lp_upolynomial_div_pseudo(&div, &rem, q, p);
        result = lp_upolynomial_is_zero(rem);
        lp_upolynomial_delete(div);
        lp_upolynomial_delete(rem);
      }
    }
  }

  if (trace_is_enabled("arithmetic")) {
    tracef("upolynomial_div_exact("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = %d\n", result);
  }

  return result;
}



void lp_upolynomial_content_Z(const lp_upolynomial_t* p, lp_integer_t* content) {

  assert(p);
  assert(p->K == lp_Z);

  unsigned i;
  lp_integer_t tmp;
  integer_construct_from_int(lp_Z, &tmp, 0);

  // GCD start
  integer_assign(lp_Z, content, &p->monomials[0].coefficient);
  if (integer_sgn(lp_Z, content) < 0) {
    integer_neg(lp_Z, &tmp, content);
    integer_swap(&tmp, content);
  }
  // GCD rest
  for (i = 1; i < p->size; ++ i) {
    integer_gcd_Z(&tmp, content, &p->monomials[i].coefficient);
    integer_swap(&tmp, content);
  }

  assert(integer_sgn(lp_Z, content) > 0);

  const lp_integer_t* lc = lp_upolynomial_lead_coeff(p);
  int sgn = integer_sgn(lp_Z, lc);
  if (sgn < 0) {
    // Content is the same sign as lc
    integer_neg(lp_Z, &tmp, content);
    integer_swap(&tmp, content);
  }

  integer_destruct(&tmp);
}

int lp_upolynomial_is_primitive(const lp_upolynomial_t* p) {
  assert(p->K == lp_Z);
  lp_integer_t content;
  integer_construct_from_int(lp_Z, &content, 0);
  lp_upolynomial_content_Z(p, &content);
  int is_primitive = integer_cmp_int(lp_Z, &content, 1) == 0 && integer_sgn(lp_Z, lp_upolynomial_lead_coeff(p)) > 0;
  integer_destruct(&content);
  return is_primitive;
}

void lp_upolynomial_make_primitive_Z(lp_upolynomial_t* p) {
  assert(p->K == lp_Z);

  lp_integer_t gcd;
  integer_construct_from_int(lp_Z, &gcd, 0);
  lp_upolynomial_content_Z(p, &gcd);

  lp_integer_t tmp;
  integer_construct_from_int(lp_Z, &tmp, 0);
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    assert(integer_divides(lp_Z, &gcd, &p->monomials[i].coefficient));
    integer_div_exact(lp_Z, &tmp, &p->monomials[i].coefficient, &gcd);
    integer_swap(&tmp, &p->monomials[i].coefficient);
  }

  integer_destruct(&gcd);
  integer_destruct(&tmp);
}

lp_upolynomial_t* lp_upolynomial_primitive_part_Z(const lp_upolynomial_t* p) {
  assert(p);
  assert(p->K == lp_Z);

  lp_upolynomial_t* result = lp_upolynomial_construct_copy(p);
  lp_upolynomial_make_primitive_Z(result);

  return result;
}

void lp_upolynomial_evaluate_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x, lp_integer_t* value) {
  const lp_int_ring_t* K = p->K;
  lp_integer_t power;
  integer_construct_from_int(lp_Z, &power, 0);

  // Compute
  integer_assign_int(lp_Z, value, 0);
  size_t i;
  for (i = 0; i < p->size; ++ i) {
    integer_pow(K, &power, x, p->monomials[i].degree);
    integer_add_mul(K, value, &p->monomials[i].coefficient, &power);
  }

  integer_destruct(&power);
}

void lp_upolynomial_evaluate_at_rational(const lp_upolynomial_t* p, const lp_rational_t* x, lp_rational_t* value) {
  assert(p->K == lp_Z);
  upolynomial_dense_t p_d;
  upolynomial_dense_construct_p(&p_d, lp_upolynomial_degree(p) + 1, p);
  upolynomial_dense_evaluate_at_rational(&p_d, x, value);
  upolynomial_dense_destruct(&p_d);
}

void lp_upolynomial_evaluate_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x, lp_dyadic_rational_t* value) {
  assert(p->K == lp_Z);
  upolynomial_dense_t p_d;
  upolynomial_dense_construct_p(&p_d, lp_upolynomial_degree(p) + 1, p);
  upolynomial_dense_evaluate_at_dyadic_rational(&p_d, x, value);
  upolynomial_dense_destruct(&p_d);
}

int lp_upolynomial_sgn_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x) {
  lp_integer_t value;
  integer_construct_from_int(p->K, &value, 0);
  lp_upolynomial_evaluate_at_integer(p, x, &value);
  int sgn = integer_sgn(p->K, &value);
  integer_destruct(&value);
  return sgn;
}

int lp_upolynomial_sgn_at_rational(const lp_upolynomial_t* p, const lp_rational_t* x) {
  lp_rational_t value;
  rational_construct(&value);
  lp_upolynomial_evaluate_at_rational(p, x, &value);
  int sgn = rational_sgn(&value);
  rational_destruct(&value);
  return sgn;
}

int lp_upolynomial_sgn_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x) {
  lp_dyadic_rational_t value;
  dyadic_rational_construct(&value);
  lp_upolynomial_evaluate_at_dyadic_rational(p, x, &value);
  int sgn = dyadic_rational_sgn(&value);
  dyadic_rational_destruct(&value);
  return sgn;
}

lp_upolynomial_t* lp_upolynomial_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q) {

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p->K == lp_Z || p->K->is_prime); // Otherwise make sure you understand what's happening

  lp_upolynomial_t* gcd = 0;

  if (lp_upolynomial_is_zero(p)) {
    gcd = lp_upolynomial_construct_copy(q);
  } else if (lp_upolynomial_is_zero(q)) {
    gcd = lp_upolynomial_construct_copy(p);
  } else if (lp_upolynomial_degree(p) < lp_upolynomial_degree(q)) {
    gcd = lp_upolynomial_gcd(q, p);
  } else {
    if (p->K == lp_Z) {
      gcd = upolynomial_gcd_heuristic(p, q, 2);
      if (!gcd) {
        gcd = upolynomial_gcd_subresultant(p, q);
      }
    } else {
      gcd = upolynomial_gcd_euclid(p, q, 0, 0);
    }
  }

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(gcd, trace_out); tracef("\n");
  }

  return gcd;
}

lp_upolynomial_t* lp_upolynomial_extended_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q, lp_upolynomial_t** u, lp_upolynomial_t** v) {

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(")\n");
  }

  assert(p->K && p->K->is_prime);
  assert(*u == 0);
  assert(*v == 0);

  lp_upolynomial_t* gcd = 0;

  if (lp_upolynomial_degree(p) < lp_upolynomial_degree(q)) {
    gcd = lp_upolynomial_extended_gcd(q, p, v, u);
  } else {
    gcd = upolynomial_gcd_euclid(p, q, u, v);
  }

  if (trace_is_enabled("gcd")) {
    tracef("upolynomial_gcd("); lp_upolynomial_print(p, trace_out); tracef(", "); lp_upolynomial_print(q, trace_out); tracef(") = "); lp_upolynomial_print(gcd, trace_out); tracef("\n");
  }

  return gcd;
}

void lp_upolynomial_solve_bezout(const lp_upolynomial_t* p, const lp_upolynomial_t* q, const lp_upolynomial_t* r, lp_upolynomial_t** u, lp_upolynomial_t** v) {
  assert(*u == 0 && *v == 0);

  lp_upolynomial_t* u1 = 0;
  lp_upolynomial_t* v1 = 0;

  // gcd = u1*p + v1*q
  lp_upolynomial_t* gcd = lp_upolynomial_extended_gcd(p, q, &u1, &v1);
  // m = r/gcd
  lp_upolynomial_t* m = lp_upolynomial_div_exact(r, gcd);
  // u = u*m
  lp_upolynomial_t* u2 = lp_upolynomial_mul(u1, m);
  // v = v*m
  lp_upolynomial_t* v2  = lp_upolynomial_mul(v1, m);
  // Fix degrees
  *u = lp_upolynomial_rem_exact(u2, q);
  *v = lp_upolynomial_rem_exact(v2, p);

  lp_upolynomial_delete(u1);
  lp_upolynomial_delete(v1);
  lp_upolynomial_delete(u2);
  lp_upolynomial_delete(v2);
  lp_upolynomial_delete(gcd);
  lp_upolynomial_delete(m);
}

lp_upolynomial_factors_t* lp_upolynomial_factor(const lp_upolynomial_t* p) {

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor("); lp_upolynomial_print(p, trace_out); tracef(")\n");
  }

  lp_upolynomial_factors_t* factors = 0;

  if (p->K == lp_Z) {
    factors = upolynomial_factor_Z(p);
  } else {
    assert(p->K->is_prime);
    factors = upolynomial_factor_Zp(p);
  }

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor("); lp_upolynomial_print(p, trace_out); tracef(") = ");
    lp_upolynomial_factors_print(factors, trace_out); tracef("\n");
  }

  return factors;
}

int lp_upolynomial_roots_count(const lp_upolynomial_t* p, const lp_rational_interval_t* ab) {
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_real_roots_count("); lp_upolynomial_print(p, trace_out); tracef(")\n");
  }
  int roots = upolynomial_roots_count_sturm(p, ab);
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_real_roots_count("); lp_upolynomial_print(p, trace_out); tracef(") => %d\n", roots);
  }
  return roots;
}

void lp_upolynomial_roots_isolate(const lp_upolynomial_t* p, lp_algebraic_number_t* roots, size_t* roots_size) {
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_roots_isolate("); lp_upolynomial_print(p, trace_out); tracef(")\n");
  }
  upolynomial_roots_isolate_sturm(p, roots, roots_size);
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_roots_isolate("); lp_upolynomial_print(p, trace_out); tracef(") => %zu\n", *roots_size);
  }
}

void lp_upolynomial_sturm_sequence(const lp_upolynomial_t* f, lp_upolynomial_t*** S, size_t* size) {
  if (trace_is_enabled("roots")) {
    tracef("upolynomial_roots_sturm_sequence("); lp_upolynomial_print(f, trace_out); tracef(")\n");
  }
  assert(f->K == lp_Z);

  size_t f_deg = lp_upolynomial_degree(f);
  upolynomial_dense_t* S_dense = (upolynomial_dense_t*) malloc((f_deg + 1)*sizeof(upolynomial_dense_t));

  upolynomial_compute_sturm_sequence(f, S_dense, size);

  (*S) = (lp_upolynomial_t**) malloc((*size)*sizeof(lp_upolynomial_t*));

  size_t i;
  for (i = 0; i < *size; ++ i) {
    (*S)[i] = upolynomial_dense_to_upolynomial(S_dense + i, lp_Z);
    upolynomial_dense_destruct(S_dense + i);
  }

  free(S_dense);
}

lp_upolynomial_t* lp_upolynomial_neg(const lp_upolynomial_t* p) {
  lp_upolynomial_t* neg = lp_upolynomial_construct_copy(p);
  lp_upolynomial_neg_in_place(neg);
  return neg;
}

void lp_upolynomial_neg_in_place(lp_upolynomial_t* p) {
  size_t i;
  for (i = 0; i < p->size; ++ i) {
    integer_neg(p->K, &p->monomials[i].coefficient, &p->monomials[i].coefficient);
  }
}

void lp_upolynomial_reverse_in_place(lp_upolynomial_t* p) {
  size_t p_deg;
  ulp_monomial_t tmp, *m_i, *m_j;

  // reverse order and update degrees
  assert(p->size > 0);
  p_deg = lp_upolynomial_degree(p);
  m_i = p->monomials;
  m_j = p->monomials + p->size - 1;
  for (; m_i <= m_j; ++ m_i, -- m_j) {
    m_i->degree = p_deg - m_i->degree;
    if (m_i < m_j) {
      m_j->degree = p_deg - m_j->degree;
      tmp = *m_i; *m_i = *m_j; *m_j = tmp;
    }
  }
}

lp_upolynomial_t* lp_upolynomial_subst_x_neg(const lp_upolynomial_t* f) {

  lp_upolynomial_t* neg = lp_upolynomial_construct_copy(f);
  size_t i;
  for (i = 0; i < neg->size; ++ i) {
    // We negate odd degrees
    if (neg->monomials[i].degree % 2) {
      integer_neg(neg->K, &neg->monomials[i].coefficient, &neg->monomials[i].coefficient);
    }
  }

  return neg;
}
