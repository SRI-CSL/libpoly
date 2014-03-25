/*
 * algebraic_number.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include "number/algebraic_number.h"
#include "polynomial_context.h"

#include "number/integer.h"
#include "interval/interval.h"
#include "interval/arithmetic.h"
#include "polynomial/coefficient.h"
#include "upolynomial/upolynomial.h"

#include "utils/debug_trace.h"

#include <assert.h>

void algebraic_number_construct(algebraic_number_t* a, upolynomial_t* f, const dyadic_interval_t* lr) {
  assert(f);
  assert(lr->a_open && lr->b_open);
  assert(upolynomial_is_primitive(f));
  a->f = f;
  dyadic_interval_construct_copy(&a->I, lr);
  a->sgn_at_a = upolynomial_sgn_at_dyadic_rational(f, &a->I.a);
  a->sgn_at_b = upolynomial_sgn_at_dyadic_rational(f, &a->I.b);
  assert(a->sgn_at_a * a->sgn_at_b < 0);
}

void algebraic_number_construct_zero(algebraic_number_t* a) {
  a->f = 0;
  dyadic_interval_construct_from_int(&a->I, 0, 0, 0, 0);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void algebraic_number_construct_copy(algebraic_number_t* a1, const algebraic_number_t* a2) {
  a1->f = a2->f ? upolynomial_construct_copy(a2->f) : 0;
  dyadic_interval_construct_copy(&a1->I, &a2->I);
  a1->sgn_at_a = a2->sgn_at_a;
  a1->sgn_at_b = a2->sgn_at_b;
}

void algebraic_number_construct_from_dyadic_rational(algebraic_number_t* a, const dyadic_rational_t* q) {
  a->f = 0;
  dyadic_interval_construct_point(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void algebraic_number_destruct(algebraic_number_t* a) {
  if (a->f) {
    upolynomial_destruct(a->f);
  }
  dyadic_interval_destruct(&a->I);
}

static inline
void algebraic_number_reduce_polynomial(const algebraic_number_t* a, const upolynomial_t* f, int sgn_at_a, int sgn_at_b) {
  assert(a->f);
  assert(a->sgn_at_a * a->sgn_at_b < 0);
  assert(sgn_at_a * sgn_at_b < 0);
  assert(upolynomial_is_primitive(f));
  algebraic_number_t* a_nonconst = (algebraic_number_t*) a;
  upolynomial_destruct(a_nonconst->f);
  a_nonconst->f = upolynomial_construct_copy(f);
}

static inline
void algebraic_number_collapse_to_point(const algebraic_number_t* a_const, const dyadic_rational_t* q) {
  assert(a_const->f);
  assert(upolynomial_sgn_at_dyadic_rational(a_const->f, q) == 0);
  // We'll modify the number so unconst it
  algebraic_number_t* a = (algebraic_number_t*) a_const;
  upolynomial_destruct(a->f);
  a->f = 0;
  dyadic_interval_collapse_to(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

static inline
void algebraic_number_refine_with_point(const algebraic_number_t* a_const, const dyadic_rational_t* q) {
  algebraic_number_t* a = (algebraic_number_t*) a_const;
  if (a->f && dyadic_interval_contains(&a->I, q)) {
    // Compute the sign at the left end point
    int a_sgn_at_q = upolynomial_sgn_at_dyadic_rational(a->f, q);
    if (a_sgn_at_q == 0) {
      algebraic_number_collapse_to_point(a, q);
      return;
    } else if (a_sgn_at_q * a->sgn_at_a > 0) {
      // We can keep I.a for a's a
      dyadic_interval_set_a(&a->I, q, 1);
    } else  {
      // We can keep I.a for a's b
      dyadic_interval_set_b(&a->I, q, 1);
    }
  }
}

/**
 * Refine the interval using one split. Returns -1 if left, +1 if right, 0 if
 * reduced to point.
 */
static inline
int algebraic_number_refine_const(const algebraic_number_t* a_const) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine(");
    algebraic_number_print(a_const, trace_out);
    tracef(")\n");
  }

  assert(a_const->f);

  int result;

  // We'll modify the number so unconst it
  algebraic_number_t* a = (algebraic_number_t*) a_const;
  // Compute the mid point
  dyadic_interval_t I_left;
  dyadic_interval_t I_right;
  dyadic_interval_construct_from_split(&I_left, &I_right, &a_const->I, 1, 1);
  // Compute the sign at the mid-point
  const dyadic_rational_t* m = &I_left.b;
  int sgn_at_m = upolynomial_sgn_at_dyadic_rational(a_const->f, m);
  if (sgn_at_m == 0) {
    // m is actually the number a1
    algebraic_number_collapse_to_point(a_const, m);
    result = 0;
  } else if (sgn_at_m * a_const->sgn_at_a > 0) {
    // a m b  or  a m b
    // + + -      - - +
    dyadic_interval_swap(&I_right, &a->I);
    result = 1;
  } else {
    // a m b  or  a m b
    // + - -      - + +
    dyadic_interval_swap(&I_left, &a->I);
    result = -1;
  }
  // Remove temp
  dyadic_interval_destruct(&I_left);
  dyadic_interval_destruct(&I_right);

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine() => ");
    algebraic_number_print(a_const, trace_out);
    tracef(", d = %d\n", result);
  }

  return result;
}

void algebraic_number_refine(algebraic_number_t* a) {
  if (a->f) {
    algebraic_number_refine_const(a);
  }
}

/**
 * The "proper" algebraic numberis always a1.
 */
int algebraic_number_cmp(const algebraic_number_t* a1, const algebraic_number_t* a2) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    algebraic_number_print(a1, trace_out);
    tracef(", ");
    algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // We only have a problem if the intervals intersect
  if (!dyadic_interval_disjunct(&a1->I, &a2->I)) {
    // First intersect the intervals. Since both intervals are open or points
    // the intersection is either open or a point
    dyadic_interval_t I;
    dyadic_interval_construct_intersection(&I, &a1->I, &a2->I);

    // Refine the interval using the intersection points
    algebraic_number_refine_with_point(a1, &I.a);
    algebraic_number_refine_with_point(a2, &I.a);
    if (!I.is_point) {
      algebraic_number_refine_with_point(a1, &I.b);
      algebraic_number_refine_with_point(a2, &I.b);
    }

    // Remove the temp intersection
    dyadic_interval_destruct(&I);
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    algebraic_number_print(a1, trace_out);
    tracef(", ");
    algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // At this point the intervals can intersect only if they are equal
  // In this case
  int equal = 0;
  if (a1->f && a2->f && dyadic_interval_equals(&a1->I, &a2->I)) {
    // If we are of the same size, and equal intervals, check for equality
    // We check by first intersecting the two
    // Both proper algebraic, intervals are equal
    // Let f = f'*gcd, g = g'*gcd
    // Gcd has a zero in the interval iff the numbers are equal and
    // we can reduce polynomials to the gcd
    upolynomial_t* gcd = upolynomial_gcd(a1->f, a2->f);
    int sgn_at_a = upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.a);
    int sgn_at_b = upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.b);
    if (sgn_at_a * sgn_at_b < 0) {
      algebraic_number_reduce_polynomial(a1, gcd, sgn_at_a, sgn_at_b);
      algebraic_number_reduce_polynomial(a2, gcd, sgn_at_a, sgn_at_b);
      equal = 1;
    } else {
      // We're not equal, so bisect away
      int d1 = 1, d2 = 1;
      while (d1 == d2 && d1 && d2) {
        // They become different when bisection goes different ways
        d1 = algebraic_number_refine_const(a1);
        d2 = algebraic_number_refine_const(a2);
      }
    }
    upolynomial_destruct(gcd);
  }

  int result;

  if (equal) {
    // We checked for equality and they are equal
    result = 0;
  } else {
    // Not equal, but disjunct, so compare the two intervals
    int cmp = dyadic_rational_ops.cmp(&a1->I.a, &a2->I.a);
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
    algebraic_number_print(a1, trace_out);
    tracef(", ");
    algebraic_number_print(a2, trace_out);
    tracef(") => %d\n", result);
  }

  return result;
}

int algebraic_number_cmp_void(const void* a1, const void* a2) {
  return algebraic_number_cmp(a1, a2);
}

int algebraic_number_print(const algebraic_number_t* a, FILE* out) {
  if (a->f == 0) {
    return dyadic_rational_ops.print(&a->I.a, out);
  } else {
    int ret = 0;
    ret += fprintf(out, "<");
    ret += upolynomial_print(a->f, out);
    ret += fprintf(out, ", ");
    ret += dyadic_interval_print(&a->I, out);
    ret += fprintf(out, ">");
    return ret;
  }
}

char* algebraic_number_to_string(const algebraic_number_t* a) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  algebraic_number_print(a, f);
  fclose(f);
  return str;
}

double algebraic_number_to_double(const algebraic_number_t* a_const) {

  // If a point, just return it's double
  if (a_const->f == 0) {
    return dyadic_rational_ops.to_double(&a_const->I.a);
  }

  // We do the necessary refinement on a copy
  algebraic_number_t a;
  algebraic_number_construct_copy(&a, a_const);

  // Refine the number until we get the desired precision
  dyadic_rational_t interval_size;
  dyadic_rational_ops.construct(&interval_size);
  dyadic_rational_ops.sub(&interval_size, &a.I.b, &a.I.a);
  if (interval_size.n < 100) {
    int iterations = 100 - interval_size.n;
    while (a.f && iterations > 0) {
      algebraic_number_refine_const(&a);
      iterations --;
    }
  }

  double result = dyadic_rational_ops.to_double(&a.I.a);

  dyadic_rational_ops.destruct(&interval_size);
  algebraic_number_destruct(&a);

  return result;
}

/** Polynomial context for resultant computation */
static const polynomial_context_t* algebraic_ctx = 0;

/** Result variable */
static variable_t var_r;

/** Variable for left operand */
static variable_t var_x;

/** Variable for right operand */
static variable_t var_y;

static
const polynomial_context_t* algebraic_pctx(void) {
  // Create the context if needed
  if (algebraic_ctx == 0) {
    // Create the variables
    variable_db_t* var_db = variable_db_ops.new();
    var_x = variable_db_ops.new_variable(var_db, "_x");
    var_y = variable_db_ops.new_variable(var_db, "_y");
    var_r = variable_db_ops.new_variable(var_db, "_r");
    // Order as Z[r, y, x]
    variable_order_simple_t* var_order = (variable_order_simple_t*) variable_order_simple_ops.variable_order_ops.new();
    variable_order_simple_ops.push(var_order, var_r);
    variable_order_simple_ops.push(var_order, var_y);
    variable_order_simple_ops.push(var_order, var_x);
    // Create the context
    algebraic_ctx = polynomial_context_ops.new(0, var_db, (variable_order_t*) var_order);
    // Detach local references
    variable_db_ops.detach(var_db);
    variable_order_simple_ops.variable_order_ops.detach((variable_order_t*) var_order);
  }
  return algebraic_ctx;
}

void filter_roots(algebraic_number_t* roots, size_t* roots_size, const dyadic_interval_t* I) {
  size_t i, to_keep;
  for (i = 0, to_keep = 0; i < *roots_size; ++ i) {
    algebraic_number_t* root = roots + i;
    if (dyadic_interval_disjunct(&root->I, I)) {
      // Remove this root
      algebraic_number_destruct(root);
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
typedef void (*construct_op_polynomial_f) (coefficient_t* op);

/** Function type called on interval operations (such as I = I1 + I2) */
typedef void (*interval_op_f) (dyadic_interval_t* I, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

/** Compute the operation */
static
void algebraic_number_op(
    algebraic_number_t* op, const algebraic_number_t* a, const algebraic_number_t* b,
    construct_op_polynomial_f construct_op,
    interval_op_f interval_op)
{
  assert(a->f && (b == 0 || b->f));

  const polynomial_context_t* ctx = algebraic_pctx();

  if (trace_is_enabled("algebraic_number")) {
    tracef("a = "); algebraic_number_print(a, trace_out); tracef("\n");
    if (b) {
      tracef("b = "); algebraic_number_print(b, trace_out); tracef("\n");
    }
    tracef("op = "); algebraic_number_print(op, trace_out); tracef("\n");
  }

  coefficient_t f_a;
  coefficient_construct_from_univariate(ctx, &f_a, a->f, var_x);

  coefficient_t f_b;
  if (b) {
    coefficient_construct_from_univariate(ctx, &f_b, b->f, var_y);
  }

  // Construct the op polynomial
  coefficient_t f_r;
  construct_op(&f_r);

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
  upolynomial_t* f = coefficient_to_univariate(ctx, &f_r);

  if (trace_is_enabled("algebraic_number")) {
    tracef("f = "); upolynomial_print(f, trace_out);
  }

  // Get the roots of f
  size_t f_roots_size = 0;
  algebraic_number_t* f_roots = malloc(sizeof(algebraic_number_t)*upolynomial_degree(f));

  upolynomial_roots_isolate(f, f_roots, &f_roots_size);

  // Interval for the result
  dyadic_interval_t I;
  dyadic_interval_construct_zero(&I);

  // Approximate the result
  while (f_roots_size > 1) {

    // Add the two intervals
    if (b) {
      interval_op(&I, &a->I, &b->I);
    } else {
      interval_op(&I, &a->I, 0);
    }

    if (trace_is_enabled("algebraic_number")) {
      tracef("a = "); algebraic_number_print(a, trace_out); tracef("\n");
      if (b) {
        tracef("b = "); algebraic_number_print(b, trace_out); tracef("\n");
      }
      tracef("I = "); dyadic_interval_print(&I, trace_out); tracef("\n");
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        tracef("f[%zu] = ", i); algebraic_number_print(f_roots + i, trace_out); tracef("\n");
      }
    }

    // Filter the roots over I
    filter_roots(f_roots, &f_roots_size, &I);

    // If more then one root, we need to refine a and b
    if (f_roots_size > 1) {
      algebraic_number_refine_const(a);
      if (b) {
        algebraic_number_refine_const(b);
      }
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        if ((f_roots + i)->f) {
          // not a point, refine it
          algebraic_number_refine_const(f_roots + i);
        }
      }
    }
  }

  assert(f_roots_size == 1);

  // Our root is the only one left
  algebraic_number_destruct(op);
  *op = *f_roots;

  if (trace_is_enabled("algebraic_number")) {
    tracef("op = "); algebraic_number_print(op, trace_out);
  }

  // Remove temps
  coefficient_destruct(&f_a);
  coefficient_destruct(&f_b);
  coefficient_destruct(&f_r);
  dyadic_interval_destruct(&I);
  free(f_roots);
}

void algebraic_number_add_construct_op(coefficient_t* f_r) {
  const polynomial_context_t* ctx = algebraic_pctx();

  // Construct the polynomial z - (x + y)
  integer_t one;
  integer_construct_from_int(Z, &one, 1);
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

void algebraic_number_add(algebraic_number_t* sum, const algebraic_number_t* a, const algebraic_number_t* b) {
  assert(a->f && b->f);
  algebraic_number_op(sum, a, b, algebraic_number_add_construct_op, dyadic_interval_add);
}

void algebraic_number_sub_construct_op(coefficient_t* f_r) {
  const polynomial_context_t* ctx = algebraic_pctx();

  // Construct the polynomial z - (x - y)
  integer_t one;
  integer_construct_from_int(Z, &one, 1);
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

void algebraic_number_sub(algebraic_number_t* sub, const algebraic_number_t* a, const algebraic_number_t* b) {
  (void) sub;
  (void) a;
  (void) b;
  assert(0);
}

void algebraic_number_mul_construct_op(coefficient_t* f_r) {
  const polynomial_context_t* ctx = algebraic_pctx();

  // Construct the polynomial z - (x*y)
  integer_t one;
  integer_construct_from_int(Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub_mul(ctx, f_r, &f_x, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

void algebraic_number_mul(algebraic_number_t* mul, const algebraic_number_t* a, const algebraic_number_t* b) {
  assert(a->f && b->f);
  algebraic_number_op(mul, a, b, algebraic_number_mul_construct_op, dyadic_interval_mul);
}

void algebraic_number_pow(algebraic_number_t* pow, const algebraic_number_t* a, unsigned n) {
  (void) pow;
  (void) a;
  (void) n;
  assert(0);

}

const algebraic_number_ops_t algebraic_number_ops = {
    algebraic_number_construct,
    algebraic_number_construct_zero,
    algebraic_number_construct_copy,
    algebraic_number_construct_from_dyadic_rational,
    algebraic_number_destruct,
    algebraic_number_cmp,
    algebraic_number_cmp_void,
    algebraic_number_print,
    algebraic_number_to_string,
    algebraic_number_to_double,
    algebraic_number_refine,
    algebraic_number_add,
    algebraic_number_sub,
    algebraic_number_mul,
    algebraic_number_pow
};

