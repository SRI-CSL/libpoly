/*
 * univariate_polynomial_buffer.c
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#include "upolynomial/upolynomial_dense.h"
#include "upolynomial/upolynomial.h"

#include "utils/debug_trace.h"

#include <malloc.h>
#include <assert.h>

static inline void upolynomial_dense_normalize(upolynomial_dense_t* p_d, int_ring K) {
  int d = p_d->size - 1;
  while (d > 0 && integer_ops.sgn(Z, p_d->coefficients + d) == 0) {
    d --;
  }
  p_d->size = d + 1;
}

static void upolynomial_dense_construct(upolynomial_dense_t* p_d, size_t capacity) {
  assert(capacity > 0);
  p_d->capacity = capacity;
  p_d->size = 1;
  p_d->coefficients = malloc(capacity*sizeof(integer_t));
  int i = 0;
  for (i = 0; i < capacity; i ++) {
    integer_ops.construct_from_int(Z, p_d->coefficients + i, 0);
  }
}

static void upolynomial_dense_construct_p(upolynomial_dense_t* p_d, size_t capacity, const upolynomial_t* p) {
  assert(capacity > upolynomial_ops.degree(p));
  upolynomial_dense_construct(p_d, capacity);
  upolynomial_ops.unpack(p, p_d->coefficients);
  p_d->size = upolynomial_ops.degree(p) + 1;
}

static void upolynomial_dense_destruct(upolynomial_dense_t* p_d) {
  int i;
  for (i = 0; i < p_d->capacity; ++ i) {
    integer_ops.destruct(p_d->coefficients + i);
  }
  free(p_d->coefficients);
}

static void upolynomial_dense_swap(upolynomial_dense_t* p_d, upolynomial_dense_t* q_d) {
  size_t tmp1;

  tmp1 = p_d->capacity;
  p_d->capacity = q_d->capacity;
  q_d->capacity = tmp1;

  tmp1 = p_d->size;
  p_d->size = q_d->size;
  q_d->size = tmp1;

  integer_t* tmp2;

  tmp2 = p_d->coefficients;
  p_d->coefficients = q_d->coefficients;
  q_d->coefficients = tmp2;
}

void upolynomial_dense_assign(upolynomial_dense_t* p_d, const upolynomial_dense_t* q_d) {
  assert(p_d->capacity >= q_d->size);
  if (p_d == q_d) return;
  size_t deg;
  for (deg = 0; deg < q_d->size; ++ deg) {
    integer_ops.assign(Z, p_d->coefficients + deg, q_d->coefficients + deg);
  }
  for (; deg < p_d->size; ++ deg) {
    integer_ops.assign_int(Z, p_d->coefficients + deg, 0);
  }
  p_d->size = q_d->size;
}

static int upolynomial_dense_is_zero(const upolynomial_dense_t* p_d) {
  return p_d->size == 1 && integer_ops.sgn(Z, p_d->coefficients) == 0;
}

const integer_t* upolynomial_dense_lead_coeff(const upolynomial_dense_t* p_d) {
  assert(p_d->size > 0);
  return p_d->coefficients + p_d->size - 1;
}

void upolynomial_dense_evaluate_at_rational(const upolynomial_dense_t* p_d, const rational_t* x, rational_t* value) {
  int i;
  rational_ops.assign_int(value, 0, 1);
  for (i = p_d->size - 1; i >= 0; -- i) {
    rational_ops.mul(value, value, x);
    rational_ops.add_integer(value, value, p_d->coefficients + i);
  }
}

void upolynomial_dense_evaluate_at_dyadic_rational(const upolynomial_dense_t* p_d, const dyadic_rational_t* x, dyadic_rational_t* value) {
  int i;
  dyadic_rational_ops.assign_int(value, 0, 0);
  for (i = p_d->size - 1; i >= 0; -- i) {
    dyadic_rational_ops.mul(value, value, x);
    dyadic_rational_ops.add_integer(value, value, p_d->coefficients + i);
  }
}

static int upolynomial_dense_sgn_at_rational(const upolynomial_dense_t* p_d, const rational_t* x) {
  rational_t value;
  rational_ops.construct(&value);
  upolynomial_dense_evaluate_at_rational(p_d, x, &value);
  int sgn = rational_ops.sgn(&value);
  rational_ops.destruct(&value);
  return sgn;
}

static int upolynomial_dense_sgn_at_dyadic_rational(const upolynomial_dense_t* p_d, const dyadic_rational_t* x) {
  dyadic_rational_t value;
  dyadic_rational_ops.construct(&value);
  upolynomial_dense_evaluate_at_dyadic_rational(p_d, x, &value);
  int sgn = dyadic_rational_ops.sgn(&value);
  dyadic_rational_ops.destruct(&value);
  return sgn;
}

static int upolynomial_dense_sgn_at_plus_inf(const upolynomial_dense_t* p_d) {
  return integer_ops.sgn(Z, p_d->coefficients + p_d->size - 1);
}

static int upolynomial_dense_sgn_at_minus_inf(const upolynomial_dense_t* p_d) {
  if (upolynomial_dense_is_zero(p_d)) {
    return 0;
  }
  if (p_d->size % 2) {
    // deg + 1 is odd, deg is even
    return integer_ops.sgn(Z, p_d->coefficients + p_d->size - 1);
  } else {
    return -integer_ops.sgn(Z, p_d->coefficients + p_d->size - 1);
  }
  return 0;
}

static void upolynomial_dense_clear(upolynomial_dense_t* p_d) {
  size_t deg;
  for (deg = 0; deg < p_d->size; ++ deg) {
    integer_ops.assign_int(Z, p_d->coefficients + deg, 0);
  }
  p_d->size = 1;
}

static upolynomial_t* upolynomial_dense_to_upolynomial(const upolynomial_dense_t* p_d, int_ring K) {
  assert(p_d->size > 0);
  upolynomial_t* result = upolynomial_ops.construct(K, p_d->size - 1, p_d->coefficients);
  return result;
}

static int upolynomial_dense_print(const upolynomial_dense_t* p_d, FILE* file) {
  int len = 0;
  int k = p_d->size - 1;
  for (; k >= 0; --k) {
    int sgn = integer_ops.sgn(Z, p_d->coefficients + k);
    if (sgn) {
      if (sgn > 0) {
        fprintf(file, "+");
      }
      len += integer_ops.print(p_d->coefficients + k, file);
      len += fprintf(file, "*x^%d ", k);
    }
  }
  return len;
}

static void upolynomial_dense_touch(upolynomial_dense_t* p_d, int_ring K, size_t degree) {
  if (degree >= p_d->size) {
    assert(degree < p_d->capacity);
    p_d->size = degree + 1;
  }
}

static void upolynomial_dense_mk_primitive_z(upolynomial_dense_t* p_d, int positive) {

  int degree = p_d->size > 0 ? p_d->size - 1 : 0;
  int lc_sgn = integer_ops.sgn(Z, p_d->coefficients + degree);

  integer_t tmp;
  integer_ops.construct_from_int(Z, &tmp, 0);

  integer_t gcd;

  // GCD start
  integer_ops.construct_copy(Z, &gcd, p_d->coefficients + degree);
  if (integer_ops.sgn(Z, &gcd) < 0) {
    integer_ops.neg(Z, &gcd, p_d->coefficients + degree);
  }
  // GCD rest
  for (-- degree; degree >= 0; -- degree) {
    if (integer_ops.sgn(Z, p_d->coefficients + degree)) {
      integer_ops.gcd_Z(&tmp, &gcd, p_d->coefficients + degree);
      integer_ops.swap(Z, &tmp, &gcd);
    }
  }
  assert(integer_ops.sgn(Z, &gcd) > 0);

  // If the lc is negative, make gcd negative
  if (positive && lc_sgn < 0) {
    integer_ops.neg(Z, &tmp, &gcd);
    integer_ops.swap(Z, &tmp, &gcd);
  }

  // Make it primitive
  if (integer_ops.cmp_int(Z, &gcd, 1)) {
    for (degree = 0; degree < p_d->size; ++ degree) {
      if (integer_ops.sgn(Z, p_d->coefficients + degree)) {
        integer_ops.div_Z(&tmp, p_d->coefficients + degree, &gcd);
        integer_ops.swap(Z, &tmp, p_d->coefficients + degree);
      }
    }
  }

  integer_ops.destruct(&tmp);
  integer_ops.destruct(&gcd);
}

static void upolynomial_dense_mult_c(upolynomial_dense_t* p_d, int_ring K, const integer_t* c) {
  assert(integer_ops.sgn(K, c));
  integer_t mult;
  integer_ops.construct_from_int(Z, &mult, 0);
  int i;
  for (i = 0; i < p_d->size; ++ i) {
    if (integer_ops.sgn(Z, p_d->coefficients + i)) {
      integer_ops.mul(K, &mult, p_d->coefficients + i, c);
      integer_ops.swap(K, &mult, p_d->coefficients + i);
    }
  }
  integer_ops.destruct(&mult);
}

static void upolynomial_dense_div_c(upolynomial_dense_t* p_d, int_ring K, const integer_t* c) {
  assert(integer_ops.sgn(K, c));
  integer_t div;
  integer_ops.construct_from_int(Z, &div, 0);
  int i;
  for (i = 0; i < p_d->size; ++ i) {
    if (integer_ops.sgn(Z, p_d->coefficients + i)) {
      integer_ops.div_exact(K, &div, p_d->coefficients + i, c);
      integer_ops.swap(K, &div, p_d->coefficients + i);
    }
  }
  integer_ops.destruct(&div);
}

static void upolynomial_dense_add_mult_p_c(upolynomial_dense_t* p_d, const upolynomial_t* p, const integer_t* c) {
  assert(integer_ops.sgn(p->K, c));
  size_t needed_degree = upolynomial_ops.degree(p);
  assert(p_d->capacity > needed_degree);
  // Add
  int i;
  for (i = 0; i < p->size; ++ i) {
    integer_ops.add_mul(p->K, &p_d->coefficients[p->monomials[i].degree], &p->monomials[i].coefficient, c);
  }
  // Adjust size
  if (needed_degree >= p_d->size) {
    p_d->size = needed_degree + 1;
  }
  upolynomial_dense_normalize(p_d, p->K);
}

static void upolynomial_dense_add_mult_p_int(upolynomial_dense_t* p_d, const upolynomial_t* p, int c) {
  assert(c);
  size_t needed_degree = upolynomial_ops.degree(p);
  assert(p_d->capacity > needed_degree);
  // Add
  int i;
  for (i = 0; i < p->size; ++ i) {
    integer_ops.add_mul_int(p->K, &p_d->coefficients[p->monomials[i].degree], &p->monomials[i].coefficient, c);
  }
  // Adjust size
  if (needed_degree >= p_d->size) {
    p_d->size = needed_degree + 1;
  }
  upolynomial_dense_normalize(p_d, p->K);
}

static void upolynomial_dense_add_mult_p_mon(upolynomial_dense_t* p_d, const upolynomial_t* p, const umonomial_t* m) {
  assert(m->degree > 0 || integer_ops.sgn(p->K, &m->coefficient));
  size_t needed_degree = upolynomial_ops.degree(p) + m->degree;
  assert(p_d->capacity > needed_degree);
  int i;
  for (i = 0; i < p->size; ++ i) {
    size_t degree = p->monomials[i].degree + m->degree;
    integer_ops.add_mul(p->K, &p_d->coefficients[degree], &p->monomials[i].coefficient, &m->coefficient);
  }
  if (needed_degree >= p_d->size) {
    p_d->size = needed_degree + 1;
  }
  upolynomial_dense_normalize(p_d, p->K);
}

static void upolynomial_dense_sub_mult_p_mon(upolynomial_dense_t* p_d, const upolynomial_t* p, const umonomial_t* m) {
  assert(m->degree > 0 || integer_ops.sgn(p->K, &m->coefficient));
  size_t needed_degree = upolynomial_ops.degree(p) + m->degree;
  int i;
  for (i = 0; i < p->size; ++ i) {
    size_t degree = p->monomials[i].degree + m->degree;
    assert(p_d->capacity > degree);
    integer_ops.sub_mul(p->K, &p_d->coefficients[degree], &p->monomials[i].coefficient, &m->coefficient);
  }
  if (needed_degree >= p_d->size) {
    p_d->size = needed_degree + 1;
  }
  upolynomial_dense_normalize(p_d, p->K);
}

static void upolynomial_dense_sub_mult_mon(upolynomial_dense_t* p_d, int_ring K, const upolynomial_dense_t* p, const umonomial_t* m) {
  assert(m->degree > 0 || integer_ops.sgn(K, &m->coefficient));

  size_t needed_size = p->size + m->degree;
  int d;
  for (d = 0; d < p->size; ++ d) {
    if (integer_ops.sgn(Z, p->coefficients + d)) {
      size_t degree = d + m->degree;
      assert(p_d->capacity > degree);
      integer_ops.sub_mul(K, &p_d->coefficients[degree], p->coefficients + d, &m->coefficient);
    }
  }

  if (needed_size > p_d->size) {
    p_d->size = needed_size;
  }
  upolynomial_dense_normalize(p_d, K);
}


void upolynomial_dense_sub_mult(upolynomial_dense_t* p_d, int_ring K, const upolynomial_dense_t* p, const upolynomial_dense_t* q) {

  if (upolynomial_dense_is_zero(p) || upolynomial_dense_is_zero(q)) {
    return;
  }

  int p_deg = p->size - 1;
  int q_deg = q->size - 1;

  assert(p_d->capacity > p_deg + q_deg);

  int i, j;
  for (i = 0; i <= p_deg; ++ i) {
    if (integer_ops.sgn(Z, p->coefficients + i)) {
      for (j = 0; j <= q_deg; ++ j) {
        if (integer_ops.sgn(Z, q->coefficients + j)) {
          integer_ops.sub_mul(K, p_d->coefficients + i + j, p->coefficients + i, q->coefficients + j);
        }
      }
    }
  }

  if (p_deg + q_deg >= p_d->size) {
    p_d->size = p_deg + q_deg + 1;
  }
  upolynomial_dense_normalize(p_d, K);
}

void upolynomial_dense_negate(upolynomial_dense_t* p_d, int_ring K) {
  int i;
  for (i = 0; i < p_d->size; ++ i) {
    integer_ops.neg(K, p_d->coefficients + i, p_d->coefficients + i);
  }
}

void upolynomial_dense_div_general(int_ring K, int exact, const upolynomial_dense_t* p, const upolynomial_dense_t* q, upolynomial_dense_t* div, upolynomial_dense_t* rem) {

  if (debug_trace_ops.is_enabled("division")) {
    tracef("upolynomial_div_general(");
    int_ring_ops.print(K, trace_out);
    tracef(", ");
    upolynomial_dense_ops.print(p, trace_out);
    tracef(", ");
    upolynomial_dense_ops.print(q, trace_out);
    tracef(")\n");
  }

  assert(p);
  assert(q);
  assert(q->size <= p->size);

  // monomial we use to multiply with
  umonomial_t m;
  umonomial_construct_from_int(Z, &m, 0, 0);

  // adjustment for the powers of div
  integer_t adjust;
  integer_ops.construct_from_int(Z, &adjust, 0);

  // Degrees of p and q
  int p_deg = p->size > 0 ? p->size - 1 : 0;
  int q_deg = q->size > 0 ? q->size - 1 : 0;

  // Copy p into rem
  upolynomial_dense_ops.assign(rem, p);
  upolynomial_dense_ops.clear(div);

  int k;
  for (k = p_deg; k >= q_deg; -- k) {
    // next coefficient to eliminate
    while (exact && k >= q_deg && integer_ops.sgn(Z, rem->coefficients + k) == 0) {
      k --;
    }
    // if found, eliminate it
    if (k >= (int) q_deg) {

      if (debug_trace_ops.is_enabled("division")) {
        tracef("q = ");
        upolynomial_dense_ops.print(q, trace_out);
        tracef("\nrem = ");
        upolynomial_dense_ops.print(rem, trace_out);
        tracef("\ndiv = ");
        upolynomial_dense_ops.print(div, trace_out);
        tracef("\n");
      }

      assert(!exact || integer_ops.divides(K, q->coefficients + q_deg, rem->coefficients + k));

      // Eliminate the coefficient using q*m
      m.degree = k - q_deg;

      // Get the quotient coefficient
      if (exact) {
        // get the quotient
        integer_ops.div_exact(K, &m.coefficient, rem->coefficients + k, q->coefficients + q_deg);
      } else {
        // rem: a*x^k, q: b*x^d, so we multiply rem with b and subtract a*q
        integer_ops.assign(Z, &m.coefficient, rem->coefficients + k);
        upolynomial_dense_ops.mult_c(rem, K, q->coefficients + q_deg);
      }

      // Do the subtraction
      if (integer_ops.sgn(Z, &m.coefficient)) {
        upolynomial_dense_ops.sub_mult_mon(rem, K, q, &m);
      }

      // Put the monomial into the division
      if (exact || !integer_ops.sgn(Z, &m.coefficient)) {
        integer_ops.swap(K, &div->coefficients[m.degree], &m.coefficient);
      } else {
        // Adjust the monomial with to be lc(q)^(k-deg_q)
        integer_ops.pow(K, &adjust, q->coefficients + q_deg, m.degree);
        integer_ops.mul(K, &div->coefficients[m.degree], &m.coefficient, &adjust);
      }
      upolynomial_dense_ops.touch(div, K, m.degree);
    }
  }

  // Free temporaries
  integer_ops.destruct(&adjust);
  umonomial_destruct(&m);

  TRACE("division", "upolynomial_div_general(): done\n");

}

void upolynomial_dense_reduce_Z(const upolynomial_dense_t* p, const upolynomial_dense_t* q, integer_t* a, upolynomial_dense_t* red) {

  integer_t lcm, red_mult;
  integer_ops.construct_from_int(Z, &lcm, 0);
  integer_ops.construct_from_int(Z, &red_mult, 0);

  if (debug_trace_ops.is_enabled("division")) {
    tracef("upolynomial_dense_reduce_Z(");
    upolynomial_dense_ops.print(p, trace_out);
    tracef(", ");
    upolynomial_dense_ops.print(q, trace_out);
    tracef(")\n");
  }

  assert(p);
  assert(q);
  assert(q->size <= p->size);

  // We start with 1
  integer_ops.assign_int(Z, a, 1);

  // monomial we use to multiply with
  umonomial_t m;
  umonomial_construct_from_int(Z, &m, 0, 0);

  // Degrees of p and q
  int p_deg = p->size > 0 ? p->size - 1 : 0;
  int q_deg = q->size > 0 ? q->size - 1 : 0;

  // Copy p into rem
  upolynomial_dense_ops.assign(red, p);

  int k;
  for (k = p_deg; k >= q_deg; -- k) {
    // next coefficient to eliminate
    while (k >= q_deg && integer_ops.sgn(Z, red->coefficients + k) == 0) {
      k --;
    }
    // if found, eliminate it
    if (k >= (int) q_deg) {

      if (debug_trace_ops.is_enabled("division")) {
        tracef("q = ");
        upolynomial_dense_ops.print(q, trace_out);
        tracef("\nred = ");
        upolynomial_dense_ops.print(red, trace_out);
        tracef("\n");
      }

      // Eliminate the coefficient using q*m
      m.degree = k - q_deg;

      // Get the quotient coefficient
      if (integer_ops.divides(Z, q->coefficients + q_deg, red->coefficients + k)) {
        integer_ops.div_Z(&m.coefficient, red->coefficients + k, q->coefficients + q_deg);
      } else{
        // red: a*x^k, q: b*x^d, so we multiply red with lcm/a and subtract
        // lcm/b multiple of q
        integer_ops.lcm_Z(&lcm, red->coefficients + k, q->coefficients + q_deg);
        TRACE("division", "lcm = %C\n", &lcm);
        integer_ops.div_exact(Z, &red_mult, &lcm, red->coefficients + k);
        TRACE("division", "red_mult = %C\n", &red_mult);
        integer_ops.div_exact(Z, &m.coefficient, &lcm, q->coefficients + q_deg);
        TRACE("division", "m.c = %C\n", &m.coefficient);
        integer_ops.mul(Z, a, a, &red_mult);
        upolynomial_dense_ops.mult_c(red, Z, &red_mult);
      }

      // Do the subtraction
      upolynomial_dense_ops.sub_mult_mon(red, Z, q, &m);

      if (debug_trace_ops.is_enabled("division")) {
        tracef("red' = ");
        upolynomial_dense_ops.print(red, trace_out);
      }
    }
  }

  // Free temporaries
  integer_ops.destruct(&lcm);
  integer_ops.destruct(&red_mult);
  umonomial_destruct(&m);

  TRACE("division", "upolynomial_dense_reduce_Z(): done\n");
}

void upolynomial_dense_derivative(int_ring K, const upolynomial_dense_t* p_d, upolynomial_dense_t* p_d_prime) {
  upolynomial_dense_clear(p_d_prime);
  int deg = p_d->size - 1;
  if (deg > 0) {
    int i;
    for (i = deg; i > 0; -- i) {
      if (integer_ops.sgn(K, &p_d->coefficients[i])) {
        integer_ops.mul_int(K, &p_d_prime->coefficients[i-1], &p_d->coefficients[i], i);
      }
    }
    upolynomial_dense_touch(p_d_prime, K, deg-1);
    upolynomial_dense_normalize(p_d_prime, K);
  }
}

const upolynomial_dense_ops_t upolynomial_dense_ops = {
    upolynomial_dense_construct,
    upolynomial_dense_construct_p,
    upolynomial_dense_destruct,
    upolynomial_dense_swap,
    upolynomial_dense_assign,
    upolynomial_dense_clear,
    upolynomial_dense_is_zero,
    upolynomial_dense_lead_coeff,
    upolynomial_dense_evaluate_at_rational,
    upolynomial_dense_evaluate_at_dyadic_rational,
    upolynomial_dense_sgn_at_rational,
    upolynomial_dense_sgn_at_dyadic_rational,
    upolynomial_dense_sgn_at_plus_inf,
    upolynomial_dense_sgn_at_minus_inf,
    upolynomial_dense_to_upolynomial,
    upolynomial_dense_print,
    upolynomial_dense_touch,
    upolynomial_dense_mk_primitive_z,
    upolynomial_dense_mult_c,
    upolynomial_dense_div_c,
    upolynomial_dense_add_mult_p_c,
    upolynomial_dense_add_mult_p_int,
    upolynomial_dense_add_mult_p_mon,
    upolynomial_dense_sub_mult_p_mon,
    upolynomial_dense_sub_mult_mon,
    upolynomial_dense_negate,
    upolynomial_dense_sub_mult,
    upolynomial_dense_div_general,
    upolynomial_dense_reduce_Z,
    upolynomial_dense_derivative
};

