/*
 * univariate_polynomial.c
 *
 *  Created on: Oct 28, 2013
 *      Author: dejan
 */

#include "upolynomial/upolynomial_internal.h"
#include "upolynomial/upolynomial_dense.h"

#include "upolynomial/gcd.h"
#include "upolynomial/factorization.h"
#include "upolynomial/root_finding.h"

#include "utils/debug_trace_internal.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

static
size_t upolynomial_degree(const upolynomial_t* p) {
  assert(p);
  assert(p->size > 0);
  return p->monomials[p->size-1].degree;
}

// Construct the polynomial, but don't construct monomials
static
upolynomial_t* upolynomial_construct_empty(int_ring K, size_t size) {
  size_t malloc_size = sizeof(upolynomial_t) + size*sizeof(umonomial_t);
  upolynomial_t* new_p = (upolynomial_t*) malloc(malloc_size);
  new_p->K = K;
  new_p->size = size;
  int_ring_ops.attach((int_ring_t*)K);
  return new_p;
}

static
upolynomial_t* upolynomial_construct(int_ring K, size_t degree, const integer_t* coefficients) {

  // Compute the needed size
  unsigned i;
  unsigned size = 0;

  for (i = 0; i <= degree; ++ i) {
    if (integer_ops.sgn(K, coefficients + i)) {
      size ++;
    }
  }

  if (size == 0) {
    size = 1;
  }

  // Allocate the memory (at least one)
  upolynomial_t* new_p = upolynomial_construct_empty(K, size > 0 ? size : 1);

  // Copy the coefficients
  size = 0;
  integer_t tmp;
  integer_ops.construct_from_int(Z, &tmp, 0);
  for (i = 0; i <= degree; ++ i) {
    integer_ops.assign(K, &tmp, coefficients + i);
    if (integer_ops.sgn(Z, &tmp)) {
      umonomial_construct(K, &new_p->monomials[size], i, &tmp);
      size ++;
    }
  }
  integer_ops.destruct(&tmp);

  // If no one constructed, make 0
  if (size == 0) {
    umonomial_construct_from_int(K, &new_p->monomials[0], 0, 0);
    size = 1;
  }

  // Return the new polynomial
  return new_p;
}

static
void upolynomial_destruct(upolynomial_t* p) {
  assert(p);
  int i = 0;
  for (i = 0; i < p->size; ++ i) {
    integer_ops.destruct(&p->monomials[i].coefficient);
  }
  int_ring_ops.detach((int_ring_t*)p->K);
  free(p);
}

static
upolynomial_t* upolynomial_construct_power(int_ring K, size_t degree, long c) {
  upolynomial_t* result = upolynomial_construct_empty(K, 1);
  integer_ops.construct_from_int(K, &result->monomials[0].coefficient, c);
  result->monomials[0].degree = degree;
  return result;
}

static
upolynomial_t* upolynomial_construct_from_int(int_ring K, size_t degree, const int* coefficients) {

  unsigned i;
  integer_t real_coefficients[degree+1];

  for (i = 0; i <= degree; ++ i) {
    integer_ops.construct_from_int(K, &real_coefficients[i], coefficients[i]);
  }

  upolynomial_t* result = upolynomial_construct(K, degree, real_coefficients);

  for (i = 0; i <= degree; ++ i) {
    integer_ops.destruct(&real_coefficients[i]);
  }

  return result;
}

static
upolynomial_t* upolynomial_construct_from_long(int_ring K, size_t degree, const long* coefficients) {

  unsigned i;
  integer_t real_coefficients[degree+1];

  for (i = 0; i <= degree; ++ i) {
    integer_ops.construct_from_int(K, &real_coefficients[i], coefficients[i]);
  }

  upolynomial_t* result = upolynomial_construct(K, degree, real_coefficients);

  for (i = 0; i <= degree; ++ i) {
    integer_ops.destruct(&real_coefficients[i]);
  }

  return result;
}

static
upolynomial_t* upolynomial_construct_copy(const upolynomial_t* p) {
  assert(p);
  // Allocate the memory
  upolynomial_t* new_p = upolynomial_construct_empty(p->K, p->size);
  // Copy the coefficients
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    umonomial_construct_copy(Z, &new_p->monomials[i], &p->monomials[i]);
  }
  // Return the new polynomial
  return new_p;
}

static
upolynomial_t* upolynomial_construct_copy_K(int_ring K, const upolynomial_t* p) {

  assert(p);
  assert(K != p->K);

  // If in Z, same, or bigger, it's just a regular copy (while swapping the rings).
  if (K == Z || (p->K != Z && integer_ops.cmp(Z, &K->M, &p->K->M) >= 0)) {
    upolynomial_t* copy = upolynomial_construct_copy(p);
    int_ring_ops.detach(copy->K);
    copy->K = K;
    int_ring_ops.attach(copy->K);
    return copy;
  }

  // Compute the needed size
  unsigned i;
  size_t size = 0;
  size_t degree = 0;

  for (i = 0; i < p->size; ++ i) {
    if (integer_ops.sgn(K, &p->monomials[i].coefficient)) {
      size ++;
      degree = p->monomials[i].degree;
    }
  }

  // At least one element
  if (degree == 0) {
    size = 1;
  }

  // Allocate the memory
  upolynomial_t* new_p = upolynomial_construct_empty(K, size);

  // Copy the coefficients
  if (degree == 0) {
    if (integer_ops.sgn(K, &p->monomials[0].coefficient)) {
      umonomial_construct_copy(K, &new_p->monomials[0], &p->monomials[0]);
    } else {
      umonomial_construct_from_int(K, &new_p->monomials[0], 0, 0);
    }
  } else {
    size = 0;
    integer_t tmp;
    integer_ops.construct_from_int(Z, &tmp, 0);
    for (i = 0; i < p->size; ++ i) {
      integer_ops.assign(K, &tmp, &p->monomials[i].coefficient);
      if (integer_ops.sgn(Z, &tmp)) {
        umonomial_construct(Z, &new_p->monomials[size], p->monomials[i].degree, &tmp);
        size ++;
      }
    }
    integer_ops.destruct(&tmp);
  }

  // Return the new polynomial
  return new_p;
}

static
void upolynomial_unpack(const upolynomial_t* p, integer_t* out) {
  assert(p);
  unsigned i; // i big step, j small step
  for (i = 0; i < p->size; ++ i) {
    integer_ops.assign(Z, &out[p->monomials[i].degree], &p->monomials[i].coefficient);
  }
}

static
int_ring upolynomial_ring(const upolynomial_t* p) {
  assert(p);
  return p->K;
}

static
void upolynomial_set_ring(upolynomial_t* p, int_ring K) {
  assert(p);
  int_ring_ops.detach(p->K);
  p->K = K;
  int_ring_ops.attach(p->K);
}

static
const integer_t* upolynomial_lead_coeff(const upolynomial_t* p) {
  assert(p);
  assert(p->size > 0);
  return &p->monomials[p->size-1].coefficient;
}

static
int upolynomial_print(const upolynomial_t* p, FILE* out) {
  assert(p);
  int len = 0;
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    if (i) {
      len += fprintf(out, " + ");
    }
    len += umonomial_print(&p->monomials[p->size-i-1], out);
  }
  len += fprintf(out, " [");
  len += int_ring_ops.print(p->K, out);
  len += fprintf(out, "]");
  return len;
}

static
char* upolynomial_to_string(const upolynomial_t* p) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  upolynomial_print(p, f);
  fclose(f);
  return str;
}

static
int upolynomial_cmp(const upolynomial_t* p, const upolynomial_t* q) {
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
    cmp = integer_ops.cmp(Z, &p->monomials[p_i].coefficient, &q->monomials[q_i].coefficient);
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

static
int upolynomial_is_zero(const upolynomial_t* p) {
  if (p->size > 1) return 0;
  if (p->monomials[0].degree > 0) return 0;
  return integer_ops.sgn(Z, &p->monomials[0].coefficient) == 0;
}

static
int upolynomial_is_one(const upolynomial_t* p) {
  if (p->size > 1) return 0;
  if (p->monomials[0].degree > 0) return 0;
  return integer_ops.cmp_int(Z, &p->monomials[0].coefficient, 1) == 0;
}

static
int upolynomial_is_monic(const upolynomial_t* p) {
  return integer_ops.cmp_int(Z, upolynomial_lead_coeff(p), 1) == 0;
}

static
upolynomial_t* upolynomial_add(const upolynomial_t* p, const upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  TRACE("arithmetic", "upolynomial_add(%P, %P)\n", p, q);

  int_ring K = p->K;

  // Get the end degree
  size_t degree = upolynomial_degree(p);
  if (upolynomial_degree(q) > degree) {
    degree = upolynomial_degree(q);
  }

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_ops.construct(&tmp, degree + 1);

  // Add the buffers
  upolynomial_dense_ops.add_mult_p_int(&tmp, p, 1);
  upolynomial_dense_ops.add_mult_p_int(&tmp, q, 1);

  // Construct the result
  upolynomial_t* result = upolynomial_dense_ops.to_upolynomial(&tmp, K);
  upolynomial_dense_ops.destruct(&tmp);

  TRACE("arithmetic", "upolynomial_add(%P, %P) = %P\n", p, q, result);

  return result;
}

static
upolynomial_t* upolynomial_sub(const upolynomial_t* p, const upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  TRACE("arithmetic", "upolynomial_sub(%P, %P)\n", p, q);

  int_ring K = p->K;

  // Get the end degree
  size_t degree = upolynomial_degree(p);
  if (upolynomial_degree(q) > degree) {
    degree = upolynomial_degree(q);
  }

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_ops.construct(&tmp, degree + 1);

  // Add the buffers
  upolynomial_dense_ops.add_mult_p_int(&tmp, p, 1);
  upolynomial_dense_ops.add_mult_p_int(&tmp, q, -1);

  // Construct the result
  upolynomial_t* result = upolynomial_dense_ops.to_upolynomial(&tmp, K);
  upolynomial_dense_ops.destruct(&tmp);

  TRACE("arithmetic", "upolynomial_sub(%P, %P) = %P\n", p, q, result);

  return result;
}

static
upolynomial_t* upolynomial_multiply_simple(const umonomial_t* m, const upolynomial_t* q) {

  assert(m);
  assert(q);

  upolynomial_t* result = upolynomial_construct_copy(q);

  int i;
  for (i = 0; i < result->size; ++ i) {
    integer_ops.mul(q->K, &result->monomials[i].coefficient, &m->coefficient, &q->monomials[i].coefficient);
    result->monomials[i].degree += m->degree;
  }

  return result;
}

static
upolynomial_t* upolynomial_multiply(const upolynomial_t* p, const upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);

  TRACE("arithmetic", "upolynomial_multiply(%P, %P)\n", p, q);

  // Take p to be the smaller size
  if (p->size > q->size) {
    return upolynomial_multiply(q, p);
  }

  // If multiplying with zero
  if (upolynomial_is_zero(p) || upolynomial_is_zero(q)) {
    return upolynomial_construct_power(p->K, 0, 0);
  }

  // Special case if possible
  if (p->K == Z && p->size == 1) {
    upolynomial_t* result = upolynomial_multiply_simple(p->monomials, q);
    TRACE("arithmetic", "upolynomial_multiply(%P, %P) = %P\n", p, q, result);
    return result;
  }

  // Max degree of the multiplications
  size_t degree = upolynomial_degree(p)+upolynomial_degree(q);

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_ops.construct(&tmp, degree + 1);

  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    upolynomial_dense_ops.add_mult_p_mon(&tmp, q, &p->monomials[i]);
  }

  // Construct the result
  upolynomial_t* result = upolynomial_dense_ops.to_upolynomial(&tmp, p->K);
  upolynomial_dense_ops.destruct(&tmp);

  TRACE("arithmetic", "upolynomial_multiply(%P, %P) = %P\n", p, q, result);

  return result;
}

static
upolynomial_t* upolynomial_multiply_c(const upolynomial_t* p, const integer_t* c) {
  assert(p);

  TRACE("arithmetic", "upolynomial_multiply_c(%P, %C)\n", p, c);

  umonomial_t m;
  umonomial_construct(p->K, &m, 0, c);
  upolynomial_t* result = upolynomial_multiply_simple(&m, p);
  umonomial_destruct(&m);

  TRACE("arithmetic", "upolynomial_multiply_c(%P, %C) = %P\n", p, c, result);

  return result;
}

static
upolynomial_t* upolynomial_pow(const upolynomial_t* p, long pow) {

  TRACE("arithmetic", "upolynomial_pow(%P, %ld)\n", p, pow);

  assert(p);
  assert(pow >= 0);

  upolynomial_t* result  = 0;

  // Anything of size 1
  if (p->size == 1) {
    result = upolynomial_construct_empty(p->K, 1);
    integer_ops.construct_from_int(Z, &result->monomials[0].coefficient, 0);
    integer_ops.pow(p->K, &result->monomials[0].coefficient, &p->monomials[0].coefficient, pow);
    result->monomials[0].degree = p->monomials[0].degree * pow;
  } else {
    // TODO: optimize
    result = upolynomial_construct_copy(p);
    int i;
    for (i = 2; i <= pow; ++ i) {
      upolynomial_t* tmp = result;
      result = upolynomial_multiply(result, p);
      upolynomial_destruct(tmp);
    }
  }

  TRACE("arithmetic", "upolynomial_pow(%P, %ld) = %P\n", p, pow, result);

  return result;
}

static
upolynomial_t* upolynomial_derivative(const upolynomial_t* p) {

  TRACE("arithmetic", "upolynomial_derivative(%P)\n", p);

  // Max degree of the derivative
  size_t degree = upolynomial_degree(p);
  if (degree > 0) degree --;

  // Ensure capacity
  upolynomial_dense_t tmp;
  upolynomial_dense_ops.construct(&tmp, degree + 1);
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    size_t degree_i = p->monomials[i].degree;
    if (degree_i > 0) {
      integer_ops.mul_int(p->K, &tmp.coefficients[degree_i-1], &p->monomials[i].coefficient, degree_i);
    }
  }
  tmp.size = degree + 1;

  // Construct the result
  upolynomial_t* result = upolynomial_dense_ops.to_upolynomial(&tmp, p->K);
  upolynomial_dense_ops.destruct(&tmp);

  TRACE("arithmetic", "upolynomial_derivative(%P) = %P\n", p, result);

  return result;
}

// Performs general division loop, result is returned in the buffers
// Buffers should be passed in unconstructed
static
void upolynomial_div_general(const upolynomial_t* p, const upolynomial_t* q, upolynomial_dense_t* div, upolynomial_dense_t* rem, int exact) {

  TRACE("division", "upolynomial_div_general(%P, %P)\n", p, q);

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(upolynomial_degree(q) <= upolynomial_degree(p));

  int_ring K = p->K;

  int p_deg = upolynomial_degree(p);
  int q_deg = upolynomial_degree(q);

  upolynomial_dense_ops.construct_p(rem, p_deg + 1, p);
  upolynomial_dense_ops.construct(div, p_deg - q_deg + 1);

  // monomial we use to multiply with
  umonomial_t m;
  umonomial_construct_from_int(Z, &m, 0, 0);

  // adjustment for the powers of div
  integer_t adjust;
  integer_ops.construct_from_int(Z, &adjust, 0);

  int k;
  for (k = p_deg; k >= q_deg; -- k) {
    // next coefficient to eliminate
    while (exact && k >= q_deg && integer_ops.sgn(Z, rem->coefficients + k) == 0) {
      k --;
    }
    // if found, eliminate it
    if (k >= q_deg) {

      if (debug_trace_ops.is_enabled("division")) {
        tracef("dividing with %P at degree %d\n", q, k);
        tracef("rem = ");
        upolynomial_dense_ops.print(rem, trace_out);
        tracef("\ndiv = ");
        upolynomial_dense_ops.print(div, trace_out);
        tracef("\n");
      }

      assert(!exact || integer_ops.divides(K, upolynomial_lead_coeff(q), rem->coefficients + k));

      // Eliminate the coefficient using q*m
      m.degree = k - q_deg;

      // Get the quotient coefficient
      if (exact) {
        // get the quotient
        integer_ops.div_exact(K, &m.coefficient, rem->coefficients + k, upolynomial_lead_coeff(q));
      } else {
        // rem: a*x^k, q: b*x^d, so we multiply rem with b and subtract a*q
        integer_ops.assign(Z, &m.coefficient, rem->coefficients + k);
        upolynomial_dense_ops.mult_c(rem, K, upolynomial_lead_coeff(q));
      }

      // Do the subtraction
      if (integer_ops.sgn(K, &m.coefficient)) {
        upolynomial_dense_ops.sub_mult_p_mon(rem, q, &m);
      }

      // Put the monomial into the division
      if (exact || !integer_ops.sgn(Z, &m.coefficient)) {
        integer_ops.swap(K, &div->coefficients[m.degree], &m.coefficient);
      } else {
        // Adjust the monomial with to be lc(q)^(k-deg_q)
        integer_ops.pow(K, &adjust, upolynomial_lead_coeff(q), m.degree);
        integer_ops.mul(K, &div->coefficients[m.degree], &m.coefficient, &adjust);
      }
      upolynomial_dense_ops.touch(div, K, m.degree);
    }
  }

  // destruct temps
  integer_ops.destruct(&adjust);
  umonomial_destruct(&m);
}

static
upolynomial_t* upolynomial_div_degrees(const upolynomial_t* p, size_t a) {

  TRACE("arithmetic", "upolynomial_div_degrees(%P, %zd)\n", p, a);

  assert(a > 1);

  upolynomial_t* result = upolynomial_construct_copy(p);
  int i;
  for (i = 0; i < result->size; ++ i) {
    assert(result->monomials[i].degree % a == 0);
    result->monomials[i].degree /= a;
  }

  TRACE("arithmetic", "upolynomial_div_degrees(%P, %zd) = %P\n", p, a, result);

  return result;
}

static
upolynomial_t* upolynomial_div_exact(const upolynomial_t* p, const upolynomial_t* q) {

  TRACE("arithmetic", "upolynomial_div_exact(%P, %P)\n", p, q);

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!upolynomial_is_zero(q));

  upolynomial_t* result = 0;

  if (upolynomial_degree(p) >= upolynomial_degree(q)) {
    int_ring K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    result = upolynomial_dense_ops.to_upolynomial(&div_buffer, K);
    upolynomial_dense_ops.destruct(&div_buffer);
    upolynomial_dense_ops.destruct(&rem_buffer);
  } else {
    // 0 polynomial
    result = upolynomial_construct_power(p->K, 0, 0);
  }

  TRACE("arithmetic", "upolynomial_div_exact(%P, %P) = %P\n", p, q, result);

  return result;
}

static
upolynomial_t* upolynomial_div_exact_c(const upolynomial_t* p, const integer_t* c) {

  TRACE("arithmetic", "upolynomial_div_exact(%P, %C)\n", p, c);

  int_ring K = p->K;

  assert(p);
  assert(integer_ops.cmp_int(K, c, 0));

  upolynomial_t* result = upolynomial_construct_empty(K, p->size);

  int i;
  for (i = 0; i < p->size; ++ i) {
    result->monomials[i].degree = p->monomials[i].degree;
    integer_ops.construct_from_int(K, &result->monomials[i].coefficient, 0);
    assert(integer_ops.divides(K, c, &p->monomials[i].coefficient));
    integer_ops.div_exact(K, &result->monomials[i].coefficient, &p->monomials[i].coefficient, c);
  }

  TRACE("arithmetic", "upolynomial_div_exact(%P, %C) = %P\n", p, c, result);

  return result;
}


static
upolynomial_t* upolynomial_rem_exact(const upolynomial_t* p, const upolynomial_t* q) {

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!upolynomial_is_zero(q));

  TRACE("arithmetic", "upolynomial_rem_exact(%P, %P)\n", p, q);

  upolynomial_t* result = 0;

  if (upolynomial_degree(p) >= upolynomial_degree(q)) {
    int_ring K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    result = upolynomial_dense_ops.to_upolynomial(&rem_buffer, K);
    upolynomial_dense_ops.destruct(&rem_buffer);
    upolynomial_dense_ops.destruct(&div_buffer);
  } else {
    // rem = p
    result = upolynomial_construct_copy(p);
  }

  TRACE("arithmetic", "upolynomial_rem_exact(%P, %P) = %P\n", p, q, result);

  return result;
}

static
void upolynomial_div_rem_exact(const upolynomial_t* p, const upolynomial_t* q,
    upolynomial_t** div, upolynomial_t** rem) {

  TRACE("arithmetic", "upolynomial_div_rem_exact(%P, %P)\n", p, q);

  assert(p);
  assert(q);
  assert(p->K == q->K);
  assert(!upolynomial_is_zero(q));
  assert(*div == 0 && *rem == 0);

  if (upolynomial_degree(p) >= upolynomial_degree(q)) {
    int_ring K = p->K;
    upolynomial_dense_t rem_buffer;
    upolynomial_dense_t div_buffer;
    upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** exact */ 1);
    *div= upolynomial_dense_ops.to_upolynomial(&div_buffer, K);
    *rem = upolynomial_dense_ops.to_upolynomial(&rem_buffer, K);
    upolynomial_dense_ops.destruct(&div_buffer);
    upolynomial_dense_ops.destruct(&rem_buffer);
  } else {
    // 0 polynomial
    *div = upolynomial_construct_power(p->K, 0, 0);
    // rem = p
    *rem = upolynomial_construct_copy(p);
  }

  TRACE("arithmetic", "upolynomial_div_exact(%P, %P) = (%P, %P)\n", p, q, div, rem);
}


static
void upolynomial_div_pseudo(upolynomial_t** div, upolynomial_t** rem, const upolynomial_t* p, const upolynomial_t* q) {

  TRACE("arithmetic", "upolynomial_div_pseudo(%P, %P)\n", p, q);

  assert(!*div);
  assert(!*rem);
  assert(p->K == q->K);
  assert(!upolynomial_is_zero(q));

  assert(upolynomial_degree(p) >= upolynomial_degree(q));

  int_ring K = p->K;

  upolynomial_dense_t rem_buffer;
  upolynomial_dense_t div_buffer;

  upolynomial_div_general(p, q, &div_buffer, &rem_buffer, /** pseudo */ 0);

  *div = upolynomial_dense_ops.to_upolynomial(&div_buffer, K);
  *rem = upolynomial_dense_ops.to_upolynomial(&rem_buffer, K);

  upolynomial_dense_ops.destruct(&div_buffer);
  upolynomial_dense_ops.destruct(&rem_buffer);

  TRACE("arithmetic", "upolynomial_div_pseudo(%P, %P) = (%P, %P)\n", p, q, *div, *rem);
}

static
int upolynomial_divides(const upolynomial_t* p, const upolynomial_t* q) {

  TRACE("arithmetic", "upolynomial_div_exact(%P, %P)\n", p, q);

  assert(p->K == q->K);

  int_ring K = p->K;

  int result = 0;

  // Check that the least power of x divides the the q
  if (p->monomials[0].degree > q->monomials[0].degree) {
    result = 0;
  } else {
    // Check that the least coefficient divides
    if (!integer_ops.divides(K, &p->monomials[0].coefficient, &q->monomials[0].coefficient)) {
      result = 0;
    } else {
      // Let's just divide
      if (K && K->is_prime) {
        upolynomial_t* rem = upolynomial_rem_exact(q, p);
        result = upolynomial_ops.is_zero(rem);
        upolynomial_ops.destruct(rem);
      } else {
        upolynomial_t* div = 0;
        upolynomial_t* rem = 0;
        upolynomial_div_pseudo(&div, &rem, q, p);
        result = upolynomial_ops.is_zero(rem);
        upolynomial_ops.destruct(div);
        upolynomial_ops.destruct(rem);
      }
    }
  }

  TRACE("arithmetic", "upolynomial_div_exact(%P, %P) = %d\n", p, q, result);

  return result;
}



static
void upolynomial_content_Z(const upolynomial_t* p, integer_t* content) {

  assert(p);
  assert(p->K == Z);

  unsigned i;
  integer_t tmp;
  integer_ops.construct_from_int(Z, &tmp, 0);

  // GCD start
  integer_ops.assign(Z, content, &p->monomials[0].coefficient);
  if (integer_ops.sgn(Z, content) < 0) {
    integer_ops.neg(Z, &tmp, content);
    integer_ops.swap(Z, &tmp, content);
  }
  // GCD rest
  for (i = 1; i < p->size; ++ i) {
    integer_ops.gcd_Z(&tmp, content, &p->monomials[i].coefficient);
    integer_ops.swap(Z, &tmp, content);
  }

  assert(integer_ops.sgn(Z, content) > 0);

  const integer_t* lc = upolynomial_lead_coeff(p);
  int sgn = integer_ops.sgn(Z, lc);
  if (sgn < 0) {
    // Content is the same sign as lc
    integer_ops.neg(Z, &tmp, content);
    integer_ops.swap(Z, &tmp, content);
  }

  integer_ops.destruct(&tmp);
}

static
int upolynomial_is_primitive(const upolynomial_t* p) {
  assert(p->K == Z);
  integer_t content;
  integer_ops.construct_from_int(Z, &content, 0);
  upolynomial_content_Z(p, &content);
  int is_primitive = integer_ops.cmp_int(Z, &content, 1) == 0 && integer_ops.sgn(Z, upolynomial_lead_coeff(p)) > 0;
  integer_ops.destruct(&content);
  return is_primitive;
}

static
void upolynomial_make_primitive_Z(upolynomial_t* p) {
  assert(p->K == Z);

  integer_t gcd;
  integer_ops.construct_from_int(Z, &gcd, 0);
  upolynomial_content_Z(p, &gcd);

  integer_t tmp;
  integer_ops.construct_from_int(Z, &tmp, 0);
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    assert(integer_ops.divides(Z, &gcd, &p->monomials[i].coefficient));
    integer_ops.div_exact(Z, &tmp, &p->monomials[i].coefficient, &gcd);
    integer_ops.swap(Z, &tmp, &p->monomials[i].coefficient);
  }

  integer_ops.destruct(&gcd);
  integer_ops.destruct(&tmp);
}

static
upolynomial_t* upolynomial_primitive_part_Z(const upolynomial_t* p) {
  assert(p);
  assert(p->K == Z);

  upolynomial_t* result = upolynomial_construct_copy(p);
  upolynomial_make_primitive_Z(result);

  return result;
}

static
void upolynomial_evaluate_at_integer(const upolynomial_t* p, const integer_t* x, integer_t* value) {
  int_ring K = p->K;
  integer_t power;
  integer_ops.construct_from_int(Z, &power, 0);

  // Compute
  integer_ops.assign_int(Z, value, 0);
  int i;
  for (i = 0; i < p->size; ++ i) {
    integer_ops.pow(K, &power, x, p->monomials[i].degree);
    integer_ops.add_mul(K, value, &p->monomials[i].coefficient, &power);
  }

  integer_ops.destruct(&power);
}

static
void upolynomial_evaluate_at_rational(const upolynomial_t* p, const rational_t* x, rational_t* value) {
  assert(p->K == Z);
  upolynomial_dense_t p_d;
  upolynomial_dense_ops.construct_p(&p_d, upolynomial_ops.degree(p) + 1, p);
  upolynomial_dense_ops.evaluate_at_rational(&p_d, x, value);
  upolynomial_dense_ops.destruct(&p_d);
}

static
void upolynomial_evaluate_at_dyadic_rational(const upolynomial_t* p, const dyadic_rational_t* x, dyadic_rational_t* value) {
  assert(p->K == Z);
  upolynomial_dense_t p_d;
  upolynomial_dense_ops.construct_p(&p_d, upolynomial_ops.degree(p) + 1, p);
  upolynomial_dense_ops.evaluate_at_dyadic_rational(&p_d, x, value);
  upolynomial_dense_ops.destruct(&p_d);
}

int upolynomial_sgn_at_integer(const upolynomial_t* p, const integer_t* x) {
  integer_t value;
  integer_ops.construct_from_int(p->K, &value, 0);
  upolynomial_evaluate_at_integer(p, x, &value);
  int sgn = integer_ops.sgn(p->K, &value);
  integer_ops.destruct(&value);
  return sgn;
}

int upolynomial_sgn_at_rational(const upolynomial_t* p, const rational_t* x) {
  rational_t value;
  rational_ops.construct(&value);
  upolynomial_evaluate_at_rational(p, x, &value);
  int sgn = rational_ops.sgn(&value);
  rational_ops.destruct(&value);
  return sgn;
}

int upolynomial_sgn_at_dyadic_rational(const upolynomial_t* p, const dyadic_rational_t* x) {
  dyadic_rational_t value;
  dyadic_rational_ops.construct(&value);
  upolynomial_evaluate_at_dyadic_rational(p, x, &value);
  int sgn = dyadic_rational_ops.sgn(&value);
  dyadic_rational_ops.destruct(&value);
  return sgn;
}

static
upolynomial_t* upolynomial_gcd(const upolynomial_t* p, const upolynomial_t* q) {

  TRACE("gcd", "upolynomial_gcd(%P, %P)\n", p, q);

  assert(p->K == Z || p->K->is_prime); // Otherwise make sure you understand what's happening

  upolynomial_t* gcd = 0;

  if (upolynomial_is_zero(p)) {
    gcd = upolynomial_construct_copy(q);
  } else if (upolynomial_is_zero(q)) {
    gcd = upolynomial_construct_copy(p);
  } else if (upolynomial_degree(p) < upolynomial_degree(q)) {
    gcd = upolynomial_gcd(q, p);
  } else {
    if (p->K == Z) {
      gcd = upolynomial_gcd_heuristic(p, q, 2);
      if (!gcd) {
        gcd = upolynomial_gcd_subresultant(p, q);
      }
    } else {
      gcd = upolynomial_gcd_euclid(p, q, 0, 0);
    }
  }

  TRACE("gcd", "upolynomial_gcd(%P, %P) = %P\n", p, q, gcd);

  return gcd;
}

upolynomial_t* upolynomial_extended_gcd(const upolynomial_t* p, const upolynomial_t* q, upolynomial_t** u, upolynomial_t** v) {

  TRACE("gcd", "upolynomial_gcd(%P, %P)\n", p, q);

  assert(p->K && p->K->is_prime);
  assert(*u == 0);
  assert(*v == 0);

  upolynomial_t* gcd = 0;

  if (upolynomial_degree(p) < upolynomial_degree(q)) {
    gcd = upolynomial_extended_gcd(q, p, v, u);
  } else {
    gcd = upolynomial_gcd_euclid(p, q, u, v);
  }

  TRACE("gcd", "upolynomial_gcd(%P, %P) = %P\n", p, q, gcd);

  return gcd;
}

void upolynomial_solve_bezout(const upolynomial_t* p, const upolynomial_t* q, const upolynomial_t* r, upolynomial_t** u, upolynomial_t** v) {
  assert(*u == 0 && *v == 0);

  upolynomial_t* u1 = 0;
  upolynomial_t* v1 = 0;

  // gcd = u1*p + v1*q
  upolynomial_t* gcd = upolynomial_extended_gcd(p, q, &u1, &v1);
  // m = r/gcd
  upolynomial_t* m = upolynomial_div_exact(r, gcd);
  // u = u*m
  upolynomial_t* u2 = upolynomial_multiply(u1, m);
  // v = v*m
  upolynomial_t* v2  = upolynomial_multiply(v1, m);
  // Fix degrees
  *u = upolynomial_rem_exact(u2, q);
  *v = upolynomial_rem_exact(v2, p);

  upolynomial_destruct(u1);
  upolynomial_destruct(v1);
  upolynomial_destruct(u2);
  upolynomial_destruct(v2);
  upolynomial_destruct(gcd);
  upolynomial_destruct(m);
}

static
upolynomial_factors_t* upolynomial_factor(const upolynomial_t* p) {

  TRACE("factorization", "upolynomial_factor(%P)\n", p);

  upolynomial_factors_t* factors = 0;

  if (p->K == Z) {
    factors = upolynomial_factor_Z(p);
  } else {
    assert(p->K->is_prime);
    factors = upolynomial_factor_Zp(p);
  }

  if (debug_trace_ops.is_enabled("factorization")) {
    tracef("upolynomial_factor(%P) = ", p);
    upolynomial_factors_ops.print(factors, trace_out);
    tracef("\n");
  }

  return factors;
}

int upolynomial_roots_count(const upolynomial_t* p, const interval_t* ab) {
  TRACE("roots", "upolynomial_real_roots_count(%P)\n", p);
  int roots = upolynomial_roots_count_sturm(p, ab);
  TRACE("roots", "upolynomial_real_roots_count(%P) => %d\n", p, roots);
  return roots;
}

void upolynomial_roots_isolate(const upolynomial_t* p, algebraic_number_t* roots, size_t* roots_size) {
  TRACE("roots", "upolynomial_roots_isolate(%P)\n", p);
  upolynomial_roots_isolate_sturm(p, roots, roots_size);
  TRACE("roots", "upolynomial_real_roots_count(%P) => %d\n", p, *roots_size);
}

void upolynomial_roots_sturm_sequence(const upolynomial_t* f, upolynomial_t*** S, size_t* size) {
  TRACE("roots", "upolynomial_roots_sturm_sequence(%P)\n", f);
  assert(f->K == Z);

  size_t f_deg = upolynomial_degree(f);
  upolynomial_dense_t* S_dense = (upolynomial_dense_t*) malloc((f_deg + 1)*sizeof(upolynomial_dense_t));

  upolynomial_compute_sturm_sequence(f, S_dense, size);

  (*S) = (upolynomial_t**) malloc((*size)*sizeof(upolynomial_t*));

  int i;
  for (i = 0; i < *size; ++ i) {
    (*S)[i] = upolynomial_dense_ops.to_upolynomial(S_dense + i, Z);
    upolynomial_dense_ops.destruct(S_dense + i);
  }

  free(S_dense);
}

const upolynomial_ops_struct upolynomial_ops = {
    upolynomial_construct,
    upolynomial_construct_power,
    upolynomial_construct_from_int,
    upolynomial_construct_from_long,
    upolynomial_construct_copy,
    upolynomial_construct_copy_K,
    upolynomial_destruct,
    upolynomial_degree,
    upolynomial_ring,
    upolynomial_set_ring,
    upolynomial_lead_coeff,
    upolynomial_unpack,
    upolynomial_print,
    upolynomial_to_string,
    upolynomial_is_zero,
    upolynomial_is_one,
    upolynomial_is_monic,
    upolynomial_is_primitive,
    upolynomial_evaluate_at_integer,
    upolynomial_evaluate_at_rational,
    upolynomial_evaluate_at_dyadic_rational,
    upolynomial_sgn_at_integer,
    upolynomial_sgn_at_rational,
    upolynomial_sgn_at_dyadic_rational,
    upolynomial_cmp,
    upolynomial_add,
    upolynomial_sub,
    upolynomial_multiply,
    upolynomial_multiply_c,
    upolynomial_pow,
    upolynomial_derivative,
    upolynomial_divides,
    upolynomial_div_degrees,
    upolynomial_div_exact,
    upolynomial_div_exact_c,
    upolynomial_rem_exact,
    upolynomial_div_rem_exact,
    upolynomial_div_pseudo,
    upolynomial_content_Z,
    upolynomial_make_primitive_Z,
    upolynomial_primitive_part_Z,
    upolynomial_gcd,
    upolynomial_extended_gcd,
    upolynomial_solve_bezout,
    upolynomial_factor,
    upolynomial_factor_square_free,
    upolynomial_roots_sturm_sequence,
    upolynomial_roots_count,
    upolynomial_roots_isolate
};

static
upolynomial_factors_t* factors_construct(void) {
  upolynomial_factors_t* f = malloc(sizeof(upolynomial_factors_t));
  integer_ops.construct_from_int(Z, &f->constant, 1);
  f->size = 0;
  f->capacity = 10;
  f->factors = calloc(f->capacity, sizeof(upolynomial_t*));
  f->multiplicities = calloc(f->capacity, sizeof(size_t));
  return f;
}

#define SWAP(type, x, y) { type tmp; tmp = x; x = y; y = tmp; }

static
void factors_swap(upolynomial_factors_t* f1, upolynomial_factors_t* f2) {
  SWAP(size_t, f1->size, f2->size);
  SWAP(size_t, f1->capacity, f2->capacity);
  SWAP(upolynomial_t**, f1->factors, f2->factors);
  SWAP(size_t*, f1->multiplicities, f2->multiplicities);
}

static
void factors_clear(upolynomial_factors_t* f) {
  int i;
  integer_ops.assign_int(Z, &f->constant, 1);
  for (i = 0; i < f->size; ++i) {
    if (f->factors[i]) {
      upolynomial_ops.destruct(f->factors[i]);
    }
  }
  f->size = 0;
}

static
void factors_destruct(upolynomial_factors_t* f, int destruct_factors) {
  if (destruct_factors) {
    factors_clear(f);
  }
  integer_ops.destruct(&f->constant);
  free(f->factors);
  free(f->multiplicities);
  free(f);
}

size_t factors_size(const upolynomial_factors_t* f) {
  return f->size;
}

upolynomial_t* factors_get_factor(upolynomial_factors_t* f, size_t i, size_t* d) {
  *d = f->multiplicities[i];
  return f->factors[i];
}

const integer_t* factors_get_constant(const upolynomial_factors_t* f) {
  return &f->constant;
}

static
void factors_add(upolynomial_factors_t* f, upolynomial_t* p, size_t d) {
  // assert(upolynomial_degree(p) > 0); (we reuse this as general sets)

  if (f->size == f->capacity) {
    f->capacity *= 2;
    f->factors = realloc(f->factors, f->capacity*sizeof(upolynomial_t*));
    f->multiplicities = realloc(f->multiplicities, f->capacity*sizeof(size_t));
  }
  f->factors[f->size] = p;
  f->multiplicities[f->size] = d;
  f->size ++;
}

int factors_print(const upolynomial_factors_t* f, FILE* out) {
  int len = 0;
  len += integer_ops.print(&f->constant, out);
  int i;
  for (i = 0; i < f->size; ++ i) {
    len += fprintf(out, " * ");
    len += fprintf(out, "[");
    len += upolynomial_ops.print(f->factors[i], out);
    len += fprintf(out, "]^%zu", f->multiplicities[i]);
  }
  return len;
}

static
int_ring factors_ring(const upolynomial_factors_t* f) {
  if (f->size == 0) {
    return Z;
  } else {
    return f->factors[0]->K;
  }
}

static
void factors_set_ring(upolynomial_factors_t* f, int_ring K) {
  int i;
  for (i = 0; i < f->size; ++ i) {
    upolynomial_ops.set_ring(f->factors[i], K);
  }
}

const upolynomial_factors_ops_t upolynomial_factors_ops = {
    factors_construct,
    factors_destruct,
    factors_clear,
    factors_swap,
    factors_size,
    factors_get_factor,
    factors_get_constant,
    factors_add,
    factors_print,
    factors_ring,
    factors_set_ring
};
