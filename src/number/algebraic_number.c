/*
 * algebraic_number.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include <interval.h>
#include <upolynomial.h>
#include <algebraic_number.h>
#include <polynomial_context.h>

#include "number/integer.h"
#include "interval/arithmetic.h"
#include "polynomial/coefficient.h"
#include "polynomial/output.h"
#include "upolynomial/output.h"

#include "utils/debug_trace.h"

#include <assert.h>

void lp_algebraic_number_construct(lp_algebraic_number_t* a, lp_upolynomial_t* f, const lp_dyadic_interval_t* lr) {
  assert(f);
  assert(lr->a_open && lr->b_open);
  assert(lp_upolynomial_is_primitive(f));
  a->f = f;
  lp_dyadic_interval_construct_copy(&a->I, lr);
  a->sgn_at_a = lp_upolynomial_sgn_at_dyadic_rational(f, &a->I.a);
  a->sgn_at_b = lp_upolynomial_sgn_at_dyadic_rational(f, &a->I.b);
  assert(a->sgn_at_a * a->sgn_at_b < 0);
}

void lp_algebraic_number_construct_zero(lp_algebraic_number_t* a) {
  a->f = 0;
  lp_dyadic_interval_construct_from_int(&a->I, 0, 0, 0, 0);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void lp_algebraic_number_construct_copy(lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2) {
  a1->f = a2->f ? lp_upolynomial_construct_copy(a2->f) : 0;
  lp_dyadic_interval_construct_copy(&a1->I, &a2->I);
  a1->sgn_at_a = a2->sgn_at_a;
  a1->sgn_at_b = a2->sgn_at_b;
}

void lp_algebraic_number_construct_from_dyadic_rational(lp_algebraic_number_t* a, const lp_dyadic_rational_t* q) {
  a->f = 0;
  lp_dyadic_interval_construct_point(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void lp_algebraic_number_destruct(lp_algebraic_number_t* a) {
  if (a->f) {
    lp_upolynomial_delete(a->f);
  }
  lp_dyadic_interval_destruct(&a->I);
}

void lp_algebraic_number_swap(lp_algebraic_number_t* a, lp_algebraic_number_t* b) {
  lp_algebraic_number_t tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline
void lp_algebraic_number_reduce_polynomial(const lp_algebraic_number_t* a, const lp_upolynomial_t* f, int sgn_at_a, int sgn_at_b) {
  __unused(sgn_at_a);
  __unused(sgn_at_b);
  assert(a->f);
  assert(a->sgn_at_a * a->sgn_at_b < 0);
  assert(sgn_at_a * sgn_at_b < 0);
  assert(lp_upolynomial_is_primitive(f));
  lp_algebraic_number_t* a_nonconst = (lp_algebraic_number_t*) a;
  lp_upolynomial_delete(a_nonconst->f);
  a_nonconst->f = lp_upolynomial_construct_copy(f);
}

static inline
void lp_algebraic_number_collapse_to_point(const lp_algebraic_number_t* a_const, const lp_dyadic_rational_t* q) {
  assert(a_const->f);
  assert(lp_upolynomial_sgn_at_dyadic_rational(a_const->f, q) == 0);
  // We'll modify the number so unconst it
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  lp_upolynomial_delete(a->f);
  a->f = 0;
  lp_dyadic_interval_collapse_to(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

static inline
void lp_algebraic_number_refine_with_point(const lp_algebraic_number_t* a_const, const lp_dyadic_rational_t* q) {
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  if (a->f && lp_dyadic_interval_contains(&a->I, q)) {
    // Compute the sign at the left end point
    int a_sgn_at_q = lp_upolynomial_sgn_at_dyadic_rational(a->f, q);
    if (a_sgn_at_q == 0) {
      lp_algebraic_number_collapse_to_point(a, q);
      return;
    } else if (a_sgn_at_q * a->sgn_at_a > 0) {
      // We can keep I.a for a's a
      lp_dyadic_interval_set_a(&a->I, q, 1);
    } else  {
      // We can keep I.a for a's b
      lp_dyadic_interval_set_b(&a->I, q, 1);
    }
  }
}

/**
 * Refine the interval using one split. Returns -1 if left, +1 if right, 0 if
 * reduced to point.
 */
static inline
int lp_algebraic_number_refine_const(const lp_algebraic_number_t* a_const) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine(");
    lp_algebraic_number_print(a_const, trace_out);
    tracef(")\n");
  }

  assert(a_const->f);

  int result;

  // We'll modify the number so unconst it
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  // Compute the mid point
  lp_dyadic_interval_t I_left;
  lp_dyadic_interval_t I_right;
  lp_dyadic_interval_construct_from_split(&I_left, &I_right, &a_const->I, 1, 1);
  // Compute the sign at the mid-point
  const lp_dyadic_rational_t* m = &I_left.b;
  int sgn_at_m = lp_upolynomial_sgn_at_dyadic_rational(a_const->f, m);
  if (sgn_at_m == 0) {
    // m is actually the number a1
    lp_algebraic_number_collapse_to_point(a_const, m);
    result = 0;
  } else if (sgn_at_m * a_const->sgn_at_a > 0) {
    // a m b  or  a m b
    // + + -      - - +
    lp_dyadic_interval_swap(&I_right, &a->I);
    result = 1;
  } else {
    // a m b  or  a m b
    // + - -      - + +
    lp_dyadic_interval_swap(&I_left, &a->I);
    result = -1;
  }
  // Remove temp
  lp_dyadic_interval_destruct(&I_left);
  lp_dyadic_interval_destruct(&I_right);

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine() => ");
    lp_algebraic_number_print(a_const, trace_out);
    tracef(", d = %d\n", result);
  }

  return result;
}

void lp_algebraic_number_refine(lp_algebraic_number_t* a) {
  if (a->f) {
    lp_algebraic_number_refine_const(a);
  }
}

/**
 * The "proper" algebraic numberis always a1.
 */
int lp_algebraic_number_cmp(const lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // We only have a problem if the intervals intersect
  if (!lp_dyadic_interval_disjunct(&a1->I, &a2->I)) {
    // First intersect the intervals. Since both intervals are open or points
    // the intersection is either open or a point
    lp_dyadic_interval_t I;
    lp_dyadic_interval_construct_intersection(&I, &a1->I, &a2->I);

    // Refine the interval using the intersection points
    lp_algebraic_number_refine_with_point(a1, &I.a);
    lp_algebraic_number_refine_with_point(a2, &I.a);
    if (!I.is_point) {
      lp_algebraic_number_refine_with_point(a1, &I.b);
      lp_algebraic_number_refine_with_point(a2, &I.b);
    }

    // Remove the temp intersection
    lp_dyadic_interval_destruct(&I);
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // At this point the intervals can intersect only if they are equal
  // In this case
  int equal = 0;
  if (a1->f && a2->f && lp_dyadic_interval_equals(&a1->I, &a2->I)) {
    // If we are of the same size, and equal intervals, check for equality
    // We check by first intersecting the two
    // Both proper algebraic, intervals are equal
    // Let f = f'*gcd, g = g'*gcd
    // Gcd has a zero in the interval iff the numbers are equal and
    // we can reduce polynomials to the gcd
    lp_upolynomial_t* gcd = lp_upolynomial_gcd(a1->f, a2->f);
    int sgn_at_a = lp_upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.a);
    int sgn_at_b = lp_upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.b);
    if (sgn_at_a * sgn_at_b < 0) {
      lp_algebraic_number_reduce_polynomial(a1, gcd, sgn_at_a, sgn_at_b);
      lp_algebraic_number_reduce_polynomial(a2, gcd, sgn_at_a, sgn_at_b);
      equal = 1;
    } else {
      // We're not equal, so bisect away
      int d1 = 1, d2 = 1;
      while (d1 == d2 && d1 && d2) {
        // They become different when bisection goes different ways
        d1 = lp_algebraic_number_refine_const(a1);
        d2 = lp_algebraic_number_refine_const(a2);
      }
    }
    lp_upolynomial_delete(gcd);
  }

  int result;

  if (equal) {
    // We checked for equality and they are equal
    result = 0;
  } else {
    // Not equal, but disjunct, so compare the two intervals
    int cmp = dyadic_rational_cmp(&a1->I.a, &a2->I.a);
    if (cmp == 0) {
      // Since they are disjunct one of them is open
      if (a1->I.a_open && !a2->I.a_open) {
        // a1  (     )
        // a2  []
        result = 1;
      } else if (!a1->I.a_open && a2->I.a_open) {
        // a1  []
        // a2  (    )
        result = -1;
      } else {
        result = cmp;
      }
    } else {
      result = cmp;
    }
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(") => %d\n", result);
  }

  return result;
}

int lp_algebraic_number_cmp_void(const void* a1, const void* a2) {
  return lp_algebraic_number_cmp(a1, a2);
}

int lp_algebraic_number_print(const lp_algebraic_number_t* a, FILE* out) {
  if (a->f == 0) {
    return dyadic_rational_print(&a->I.a, out);
  } else {
    int ret = 0;
    ret += fprintf(out, "<");
    ret += lp_upolynomial_print(a->f, out);
    ret += fprintf(out, ", ");
    ret += lp_dyadic_interval_print(&a->I, out);
    ret += fprintf(out, ">");
    return ret;
  }
}

char* lp_algebraic_number_to_string(const lp_algebraic_number_t* a) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_algebraic_number_print(a, f);
  fclose(f);
  return str;
}

double lp_algebraic_number_to_double(const lp_algebraic_number_t* a_const) {

  // If a point, just return it's double
  if (a_const->f == 0) {
    return dyadic_rational_to_double(&a_const->I.a);
  }

  // We do the necessary refinement on a copy
  lp_algebraic_number_t a;
  lp_algebraic_number_construct_copy(&a, a_const);

  // Refine the number until we get the desired precision
  lp_dyadic_rational_t interval_size;
  dyadic_rational_construct(&interval_size);
  dyadic_rational_sub(&interval_size, &a.I.b, &a.I.a);
  if (interval_size.n < 100) {
    int iterations = 100 - interval_size.n;
    while (a.f && iterations > 0) {
      lp_algebraic_number_refine_const(&a);
      iterations --;
    }
  }

  double result = dyadic_rational_to_double(&a.I.a);

  dyadic_rational_destruct(&interval_size);
  lp_algebraic_number_destruct(&a);

  return result;
}

/** Polynomial context for resultant computation */
static const lp_polynomial_context_t* algebraic_ctx = 0;

/** Result variable */
static lp_variable_t var_r;

/** Variable for left operand */
static lp_variable_t var_x;

/** Variable for right operand */
static lp_variable_t var_y;

static
const lp_polynomial_context_t* lp_algebraic_pctx(void) {
  // Create the context if needed
  if (algebraic_ctx == 0) {
    // Create the variables
    lp_variable_db_t* var_db = lp_variable_db_new();
    var_x = lp_variable_db_new_variable(var_db, "_x");
    var_y = lp_variable_db_new_variable(var_db, "_y");
    var_r = lp_variable_db_new_variable(var_db, "_r");
    // Order as Z[r, y, x]
    lp_variable_order_t* var_order = lp_variable_order_new();
    lp_variable_order_push(var_order, var_r);
    lp_variable_order_push(var_order, var_y);
    lp_variable_order_push(var_order, var_x);
    // Create the context
    algebraic_ctx = lp_polynomial_context_new(0, var_db, var_order);
    // Detach local references
    lp_variable_db_detach(var_db);
    lp_variable_order_detach(var_order);
  }
  return algebraic_ctx;
}

void filter_roots(lp_algebraic_number_t* roots, size_t* roots_size, const lp_dyadic_interval_t* I) {
  size_t i, to_keep;
  for (i = 0, to_keep = 0; i < *roots_size; ++ i) {
    lp_algebraic_number_t* root = roots + i;
    if (lp_dyadic_interval_disjunct(&root->I, I)) {
      // Remove this root
      lp_algebraic_number_destruct(root);
    } else {
      // Keep it
      if (i > to_keep) {
        *(roots + to_keep) = *root;
      }
      to_keep ++;
    }
  }
  *roots_size = to_keep;
}

/** Function type called on coefficient traversal (such as r - (x + y)) */
typedef void (*construct_op_polynomial_f) (coefficient_t* op, void* data);

/** Function type called on interval operations (such as I = I1 + I2) */
typedef void (*interval_op_f) (lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data);

/** Compute the operation */
static
void lp_algebraic_number_op(
    lp_algebraic_number_t* op, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b,
    construct_op_polynomial_f construct_op,
    interval_op_f interval_op,
    void* data)
{
  assert(a->f && (b == 0 || b->f));

  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  if (trace_is_enabled("algebraic_number")) {
    tracef("a = "); lp_algebraic_number_print(a, trace_out); tracef("\n");
    if (b) {
      tracef("b = "); lp_algebraic_number_print(b, trace_out); tracef("\n");
    }
    tracef("op = "); lp_algebraic_number_print(op, trace_out); tracef("\n");
  }

  coefficient_t f_a;
  coefficient_construct_from_univariate(ctx, &f_a, a->f, var_x);

  coefficient_t f_b;
  if (b) {
    coefficient_construct_from_univariate(ctx, &f_b, b->f, var_y);
  }

  // Construct the op polynomial
  coefficient_t f_r;
  construct_op(&f_r, data);

  if (trace_is_enabled("algebraic_number")) {
    tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
  }

  // Compute the resultant
  coefficient_resultant(ctx, &f_r, &f_r, &f_a);
  if (trace_is_enabled("algebraic_number")) {
    tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
  }
  if (b) {
    coefficient_resultant(ctx, &f_r, &f_r, &f_b);
    if (trace_is_enabled("algebraic_number")) {
      tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
    }
  }

  // Resultant polynomial captures the result
  lp_upolynomial_t* f = coefficient_to_univariate(ctx, &f_r);

  if (trace_is_enabled("algebraic_number")) {
    tracef("f = "); lp_upolynomial_print(f, trace_out);
  }

  // Get the roots of f
  size_t f_roots_size = 0;
  lp_algebraic_number_t* f_roots = malloc(sizeof(lp_algebraic_number_t)*lp_upolynomial_degree(f));
  lp_upolynomial_roots_isolate(f, f_roots, &f_roots_size);
  lp_upolynomial_delete(f);

  // Interval for the result
  lp_dyadic_interval_t I;
  lp_dyadic_interval_construct_zero(&I);

  // Approximate the result
  while (f_roots_size > 1) {

    // Add the two intervals
    if (b) {
      interval_op(&I, &a->I, &b->I, data);
    } else {
      interval_op(&I, &a->I, 0, data);
    }

    if (trace_is_enabled("algebraic_number")) {
      tracef("a = "); lp_algebraic_number_print(a, trace_out); tracef("\n");
      if (b) {
        tracef("b = "); lp_algebraic_number_print(b, trace_out); tracef("\n");
      }
      tracef("I = "); lp_dyadic_interval_print(&I, trace_out); tracef("\n");
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        tracef("f[%zu] = ", i); lp_algebraic_number_print(f_roots + i, trace_out); tracef("\n");
      }
    }

    // Filter the roots over I
    filter_roots(f_roots, &f_roots_size, &I);

    // If more then one root, we need to refine a and b
    if (f_roots_size > 1) {
      lp_algebraic_number_refine_const(a);
      if (b) {
        lp_algebraic_number_refine_const(b);
      }
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        if ((f_roots + i)->f) {
          // not a point, refine it
          lp_algebraic_number_refine_const(f_roots + i);
        }
      }
    }
  }

  assert(f_roots_size == 1);

  // Our root is the only one left
  lp_algebraic_number_destruct(op);
  *op = *f_roots;

  if (trace_is_enabled("algebraic_number")) {
    tracef("op = "); lp_algebraic_number_print(op, trace_out);
  }

  // Remove temps
  coefficient_destruct(&f_a);
  if (b) {
    coefficient_destruct(&f_b);
  }
  coefficient_destruct(&f_r);
  lp_dyadic_interval_destruct(&I);
  free(f_roots);
}

static
void lp_algebraic_number_add_construct_op(coefficient_t* f_r, void* data) {
  __unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x + y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  coefficient_sub(ctx, f_r, f_r, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

static
void lp_algebraic_number_add_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __unused(data);
  dyadic_interval_add(I, I1, I2);
}

void lp_algebraic_number_add(lp_algebraic_number_t* sum, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  assert(a->f && b->f);
  lp_algebraic_number_op(sum, a, b, lp_algebraic_number_add_construct_op, lp_algebraic_number_add_interval_op, 0);
}

void lp_algebraic_number_sub_construct_op(coefficient_t* f_r, void* data) {
  __unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x - y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  coefficient_add(ctx, f_r, f_r, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

void lp_algebraic_number_sub_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __unused(data);
  dyadic_interval_sub(I, I1, I2);
}

void lp_algebraic_number_sub(lp_algebraic_number_t* sub, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  assert(a->f && b->f);
  lp_algebraic_number_op(sub, a, b, lp_algebraic_number_sub_construct_op, lp_algebraic_number_sub_interval_op, 0);
}

void lp_algebraic_number_neg(lp_algebraic_number_t* neg, const lp_algebraic_number_t* a) {
  assert(a->f);
  lp_upolynomial_t* f_neg_x = lp_upolynomial_subst_x_neg(a->f);
  if (integer_sgn(lp_Z, lp_upolynomial_lead_coeff(f_neg_x)) < 0) {
    lp_upolynomial_neg_in_place(f_neg_x);
  }
  lp_dyadic_interval_t I_neg; // To destroy
  lp_dyadic_interval_construct_copy(&I_neg, &a->I);
  dyadic_interval_neg(&I_neg, &I_neg);

  lp_algebraic_number_t result; // To destroy
  lp_algebraic_number_construct(&result, f_neg_x, &I_neg);
  lp_algebraic_number_swap(&result, neg);

  lp_algebraic_number_destruct(&result);
  lp_dyadic_interval_destruct(&I_neg);
}

void lp_algebraic_number_mul_construct_op(coefficient_t* f_r, void* data) {
  __unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x*y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub_mul(ctx, f_r, &f_x, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

void lp_algebraic_number_mul_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __unused(data);
  dyadic_interval_mul(I, I1, I2);
}

void lp_algebraic_number_mul(lp_algebraic_number_t* mul, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  assert(a->f && b->f);
  lp_algebraic_number_op(mul, a, b, lp_algebraic_number_mul_construct_op, lp_algebraic_number_mul_interval_op, 0);
}

void lp_algebraic_number_pow_construct_op(coefficient_t* f_r, void* data) {
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  unsigned n = *((unsigned*) data);

  // Construct the polynomial z - x**n
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, n);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
}

void lp_algebraic_number_pow_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  assert(I2 == 0);
  __unused(I2);
  unsigned n = *((unsigned*) data);
  dyadic_interval_pow(I, I1, n);
}

void lp_algebraic_number_pow(lp_algebraic_number_t* pow, const lp_algebraic_number_t* a, unsigned n) {
  lp_algebraic_number_op(pow, a, 0, lp_algebraic_number_pow_construct_op, lp_algebraic_number_pow_interval_op, &n);
}
