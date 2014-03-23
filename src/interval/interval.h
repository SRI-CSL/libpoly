/*
 * interval_internal.h
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#pragma once

#include <interval.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

#include <assert.h>
#include <limits.h>

static inline
void interval_construct(interval_t* I,
    const rational_t* a, int a_open,
    const rational_t* b, int b_open)
{
  int cmp = rational_cmp(a, b);
  assert(cmp <= 0);
  rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    rational_construct_copy(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    I->is_point = 1;
    assert(!a_open && !b_open);
    I->a_open = I->b_open = 0;
  }
}

static inline
void interval_construct_point(interval_t* I, const rational_t* a)
{
  rational_construct_copy(&I->a, a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

static inline
void interval_construct_zero(interval_t* I) {
  rational_construct(&I->a);
  I->a_open = I->b_open = 0;
  I->is_point = 1;
}

static inline
void interval_construct_copy(interval_t* I, const interval_t* from) {
  rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

static inline
void interval_construct_from_dyadic(interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  rational_construct_from_dyadic(&I->a, a);
  if (cmp != 0) {
    rational_construct_from_dyadic(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

static inline
void interval_construct_from_int(interval_t* I,
    long a, int a_open,
    long b, int b_open)
{
  assert(a <= b);
  rational_construct_from_int(&I->a, a, 0);
  if (a != b) {
    rational_construct_from_int(&I->b, b, 0);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

static inline
void interval_construct_from_integer(interval_t* I,
    const integer_t* a, int a_open,
    const integer_t* b, int b_open)
{
  int cmp = integer_cmp(Z, a, b);
  assert(cmp <= 0);
  rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    rational_construct_from_integer(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

static inline
void interval_destruct(interval_t* I) {
  rational_destruct(&I->a);
  if (!I->is_point) {
    rational_destruct(&I->b);
  }
}

static inline
void interval_swap(interval_t* I1, interval_t* I2) {
  interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

/**
 * Returns the sign of the interval. The sign is 0 if the interval contains
 * 0, otherwise it's the sign of all of it's points.
 */
static inline
int interval_sgn(const interval_t* I) {
  int a_sgn = rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = rational_sgn(&I->b);

  if (a_sgn < 0 && b_sgn > 0) {
    // Definitively contains 0
    return 0;
  }
  if (a_sgn == 0) {
    if (!I->a_open) {
      // Left contains zero
      return 0;
    } else {
      return 1;
    }
  }
  if (b_sgn == 0) {
    if (!I->b_open) {
      // Right contains zero
      return 0;
    } else {
      return -1;
    }
  }

  // Doesn't contain zero
  if (a_sgn < 0) {
    return -1;
  }

  assert(b_sgn > 0);
  return 1;
}

static inline
int interval_endpoint_lt(const rational_t* a, int a_open, const rational_t*b, int b_open) {
  int cmp = rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

static inline
void interval_mul(interval_t* P, const interval_t* I1, const interval_t* I2) {
  if (I1->is_point) {
    if (I2->is_point) {
      // Just multiply the points
      rational_mul(&P->a, &I1->a, &I2->a);
      if (!P->is_point) {
        rational_destruct(&P->b);
        P->is_point = 1;
      }
      P->a_open = P->b_open = 0;
    } else {
      // Depending on the sign of a, we might have to flip
      int a_sgn = rational_sgn(&I1->a);
      if (a_sgn == 0) {
        // It's just 0
        if (!P->is_point) {
          rational_destruct(&P->b);
          P->is_point = 1;
        }
        P->a_open = P->b_open = 0;
        rational_assign_int(&P->a, 0, 1);
      } else if (a_sgn > 0) {
        // Regular multiplication
        if (P->is_point) {
          rational_construct(&P->b);
          P->is_point = 0;
        }
        P->a_open = I2->a_open;
        P->b_open = I2->b_open;
        rational_mul(&P->a, &I1->a, &I2->a);
        rational_mul(&P->b, &I1->a, &I2->b);
      } else {
        // Multiplying with a negative, flip the edges
        if (P->is_point) {
          rational_construct(&P->b);
          P->is_point = 0;
        }
        P->a_open = I2->b_open;
        P->b_open = I2->a_open;
        rational_mul(&P->a, &I1->a, &I2->b);
        rational_mul(&P->b, &I1->a, &I2->a);
      }
    }
  } else if (I2->is_point) {
    interval_mul(P, I2, I1);
  } else {
    if (P->is_point) {
      rational_construct(&P->b);
      P->is_point = 0;
    }

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //         = { x*y |

    interval_t result;
    interval_construct_zero(&result);

    rational_t tmp;
    rational_construct(&tmp);

    // I1.a x I2.a
    rational_mul(&result.a, &I1->a, &I2->a);
    rational_construct_copy(&result.b, &result.a);
    result.a_open = result.b_open = I1->a_open || I2->a_open;
    result.is_point = 0;

    // I1.a x I2.b
    int tmp_open = I1->a_open || I2->b_open;
    rational_mul(&tmp, &I1->a, &I2->b);
    if (interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    tmp_open = I1->b_open || I2->a_open;
    rational_mul(&tmp, &I1->b, &I2->a);
    if (interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    tmp_open = I1->b_open || I2->b_open;
    rational_mul(&tmp, &I1->b, &I2->b);
    if (interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    interval_swap(&result, P);
    interval_destruct(&result);
    rational_destruct(&tmp);
  }
}

static inline
void interval_add(interval_t* S, const interval_t* I1, const interval_t* I2) {
  // [a, b] + [c, d] = [a + c, b + d]
  interval_t result;
  rational_construct(&result.a);
  rational_construct(&result.b);
  rational_add(&result.a, &I1->a, &I2->a);
  rational_add(&result.b, &I1->b, &I2->b);
  result.a_open = I1->a_open || I2->a_open;
  result.b_open = I1->b_open || I2->b_open;
  interval_swap(&result, S);
  interval_destruct(&result);
}

static inline
void interval_pow(interval_t* P, const interval_t* I, unsigned n) {
  if (n == 0) {
    // I^0 = [1]
    if (!P->is_point) {
      P->is_point = 1;
      rational_destruct(&P->b);
    }
    rational_assign_int(&P->a, 1, 1);
    P->a_open = 0;
    P->b_open = 0;
  } else if (I->is_point) {
    // Plain power
    if (!P->is_point) {
      rational_destruct(&P->b);
      P->is_point = 1;
      P->a_open = P->b_open = 0;
    }
    rational_pow(&P->a, &I->a, n);
  } else {
    if (P->is_point) {
      P->is_point = 0;
      rational_construct(&P->b);
    }
    if (n % 2) {
      // For odd powers we are monotonic, i.e. [a, b]^n = [a^n, b^n]
      P->a_open = I->a_open;
      P->b_open = I->b_open;
      rational_pow(&P->a, &I->a, n);
      rational_pow(&P->b, &I->b, n);
    } else {
      // Even powers depend on whether 0 is in the interval
      int sgn = interval_sgn(I);
      rational_pow(&P->a, &I->a, n);
      rational_pow(&P->b, &I->b, n);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        if (interval_endpoint_lt(&P->b, I->b_open, &P->a, I->a_open)) {
          rational_swap(&P->b, &P->a);
          P->b_open = I->a_open;
        }
        rational_assign_int(&P->a, 0, 1);
        P->a_open = 0;
      } else if (sgn > 0) {
        // P = I^n
        P->a_open = I->a_open;
        P->b_open = I->b_open;
      } else {
        // negative turns positive, so we flip
        rational_swap(&P->a, &P->b);
        P->a_open = I->b_open;
        P->b_open = I->a_open;
      }
    }
  }
}

int interval_print(const interval_t* I, FILE* out);

static inline
void dyadic_interval_construct(dyadic_interval_t* I,
    const dyadic_rational_t* a, int a_open,
    const dyadic_rational_t* b, int b_open)
{
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_copy(&I->b, b);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

static inline
void dyadic_interval_construct_zero(dyadic_interval_t* I) {
  dyadic_rational_construct(&I->a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

static inline
void dyadic_interval_construct_copy(dyadic_interval_t* I, const dyadic_interval_t* from) {
  dyadic_rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    dyadic_rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

static inline
void dyadic_interval_construct_from_int(dyadic_interval_t* I,
    long a, int a_open,
    long b, int b_open)
{
  assert(a <= b);
  dyadic_rational_construct_from_int(&I->a, a, 0);
  if (a != b) {
    dyadic_rational_construct_from_int(&I->b, b, 0);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }

}

static inline
void dyadic_interval_construct_from_integer(dyadic_interval_t* I,
    const integer_t* a, int a_open,
    const integer_t* b, int b_open)
{
  int cmp = integer_cmp(Z, a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_from_integer(&I->b, b);
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

static inline
void dyadic_interval_construct_point(dyadic_interval_t* I, const dyadic_rational_t* q) {
  dyadic_rational_construct_copy(&I->a, q);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

static inline
void dyadic_interval_construct_from_split(dyadic_interval_t* I_left, dyadic_interval_t* I_right, const dyadic_interval_t* I, int left_open, int right_open) {
  assert(!I->is_point);
  dyadic_rational_t m;
  dyadic_rational_construct(&m);
  dyadic_rational_add(&m, &I->a, &I->b);
  dyadic_rational_div_2exp(&m, &m, 1);
  dyadic_interval_construct(I_left, &I->a, I->a_open, &m, left_open);
  dyadic_interval_construct(I_right, &m, right_open, &I->b, I->b_open);
  dyadic_rational_destruct(&m);
}

static inline
int dyadic_interval_contains(const dyadic_interval_t* I, const dyadic_rational_t* q) {
  int cmp_a = dyadic_rational_cmp(&I->a, q);
  if (I->is_point) {
    return cmp_a == 0;
  }
  if (I->a_open && !(cmp_a < 0)) return 0;
  if (!I->a_open && !(cmp_a <= 0)) return 0;
  int cmp_b = dyadic_rational_cmp(q, &I->b);
  if (I->b_open && !(cmp_b < 0)) return 0;
  if (!I->b_open && !(cmp_b <= 0)) return 0;
  return 1;
}

static inline
void dyadic_interval_construct_intersection(dyadic_interval_t* I, const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point) {
    assert(dyadic_interval_contains(I2, &I1->a));
    dyadic_interval_construct_copy(I, I1);
  } else if (I2->is_point) {
    assert(dyadic_interval_contains(I1, &I2->a));
    dyadic_interval_construct_copy(I, I2);
  } else {
    // (   [  )   ]
    int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
    const dyadic_rational_t* max_a = cmp_a < 0 ? &I2->a : &I1->a;
    int a_open;
    if (cmp_a == 0) {
      a_open = I1->a_open || I2->a_open;
    } else if (cmp_a < 0) {
      a_open = I2->a_open;
    } else {
      a_open = I1->a_open;
    }

    int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
    const dyadic_rational_t* min_b = cmp_b < 0 ? &I1->b : &I2->b;
    int b_open;
    if (cmp_b == 0) {
      b_open = I1->b_open || I2->b_open;
    } else if (cmp_b < 0) {
      b_open = I1->b_open;
    } else {
      b_open = I2->b_open;
    }

    dyadic_interval_construct(I, max_a, a_open, min_b, b_open);
  }
}

static inline
void dyadic_interval_destruct(dyadic_interval_t* I) {
  dyadic_rational_destruct(&I->a);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
}

static inline
void dyadic_interval_assign(dyadic_interval_t* I, const dyadic_interval_t* from) {
  if (I != from) {
    if (I->is_point) {
      if (from->is_point) {
        dyadic_rational_assign(&I->a, &from->a);
      } else {
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_construct_copy(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    } else {
      if (from->is_point) {
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_destruct(&I->b);
        I->a_open = I->b_open = 0;
        I->is_point = 1;
      } else {
        // both intervals
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_assign(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    }
  }
}

static inline
void dyadic_interval_swap(dyadic_interval_t* I1, dyadic_interval_t* I2) {
  dyadic_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

static inline
int dyadic_interval_sgn(const dyadic_interval_t* I) {
  int a_sgn = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = dyadic_rational_sgn(&I->b);
  if (a_sgn < 0 && b_sgn <= 0) {
    return -1;
  } else if (b_sgn > 0 && a_sgn >= 0) {
    return 1;
  } else {
    // a_sgn <= 0, b_sgn >= 0, contains 0
    return 0;
  }
}

static inline
void dyadic_interval_collapse_to(dyadic_interval_t* I, const dyadic_rational_t* q) {
  dyadic_rational_assign(&I->a, q);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

static inline
void dyadic_interval_set_a(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(a, &I->b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->a, a);
    I->a_open = a_open;
  } else {
    assert(!a_open && !I->b_open);
    dyadic_interval_collapse_to(I, a);
  }
}

static inline
void dyadic_interval_set_b(dyadic_interval_t* I, const dyadic_rational_t* b, int b_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(&I->a, b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->b, b);
    I->b_open = b_open;
  } else {
    assert(!I->a_open && !b_open);
    dyadic_interval_collapse_to(I, b);
  }
}

static inline
int dyadic_interval_equals(const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point && !I2->is_point) {
    return 0;
  }
  if (!I1->is_point && I2->is_point) {
    return 0;
  }
  int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
  if (I1->is_point) {
    assert(I2->is_point);
    return cmp_a == 0;
  }
  if (cmp_a != 0 || ((!I1->a_open) != (!I2->a_open))) return 0;
  int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
  if (cmp_b != 0 || ((!I1->b_open) != (!I2->b_open))) return 0;
  return 1;
}

static inline
int dyadic_interval_disjunct(const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point) {
    return !dyadic_interval_contains(I2, &I1->a);
  }
  if (I2->is_point) {
    return !dyadic_interval_contains(I1, &I2->a);
  }
  // ( I1 ) ( I2 ) ?
  int cmp1 = dyadic_rational_cmp(&I1->b, &I2->a);
  if (cmp1 < 0) return 1;
  else if (cmp1 == 0 && (I1->b_open || I2->a_open)) return 1;
  // ( I2 ) ( I1 ) ?
  int cmp2 = dyadic_rational_cmp(&I2->b, &I1->a);
  if (cmp2 < 0) return 1;
  else if (cmp2 == 0 && (I2->b_open || I1->a_open)) return 1;
  return 0;
}

static inline
void dyadic_interval_scale(dyadic_interval_t* I, int n) {
  assert(!I->is_point);
  if (n > 0) {
    dyadic_rational_mul_2exp(&I->a, &I->a, n);
    dyadic_rational_mul_2exp(&I->b, &I->b, n);
  } else {
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
  }
}

static inline
void dyadic_interval_add(dyadic_interval_t* S, const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  // [a, b] + [c, d] = [a + c, b + d]
  if (I1->is_point && I2->is_point) {
    if (!S->is_point) {
      dyadic_rational_destruct(&S->b);
    }
    dyadic_rational_add(&S->a, &I1->a, &I2->a);
    S->b_open = S->a_open = 0;
    S->is_point = 1;
    return;
  }

  if (I2->is_point) {
    // Reuse symmetry
    dyadic_interval_add(S, I2, I1);
    return;
  }

  if (I1->is_point) {
    // Just shift by I1->a
    dyadic_interval_assign(S, I2);
    dyadic_rational_add(&S->a, &S->a, &I1->a);
    dyadic_rational_add(&S->b, &S->b, &I1->a);
    return;
  }

  // Both non-points
  dyadic_interval_t result;
  dyadic_rational_construct(&result.a);
  dyadic_rational_construct(&result.b);
  dyadic_rational_add(&result.a, &I1->a, &I2->a);
  dyadic_rational_add(&result.b, &I1->b, &I2->b);
  result.a_open = I1->a_open || I2->a_open;
  result.b_open = I1->b_open || I2->b_open;
  result.is_point = 0;
  dyadic_interval_swap(&result, S);
  dyadic_interval_destruct(&result);
}

static inline
int dyadic_interval_endpoint_lt(const dyadic_rational_t* a, int a_open, const dyadic_rational_t*b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

static inline
void dyadic_interval_mul(dyadic_interval_t* P, const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point) {
    if (I2->is_point) {
      // Just multiply the points
      dyadic_rational_mul(&P->a, &I1->a, &I2->a);
      if (!P->is_point) {
        dyadic_rational_destruct(&P->b);
        P->is_point = 1;
      }
      P->a_open = P->b_open = 0;
    } else {
      // Depending on the sign of a, we might have to flip
      int a_sgn = dyadic_rational_sgn(&I1->a);
      if (a_sgn == 0) {
        // It's just 0
        if (!P->is_point) {
          dyadic_rational_destruct(&P->b);
          P->is_point = 1;
        }
        P->a_open = P->b_open = 0;
        dyadic_rational_assign_int(&P->a, 0, 1);
      } else if (a_sgn > 0) {
        // Regular multiplication
        if (P->is_point) {
          dyadic_rational_construct(&P->b);
          P->is_point = 0;
        }
        P->a_open = I2->a_open;
        P->b_open = I2->b_open;
        dyadic_rational_mul(&P->a, &I1->a, &I2->a);
        dyadic_rational_mul(&P->b, &I1->a, &I2->b);
      } else {
        // Multiplying with a negative, flip the edges
        if (P->is_point) {
          dyadic_rational_construct(&P->b);
          P->is_point = 0;
        }
        P->a_open = I2->b_open;
        P->b_open = I2->a_open;
        dyadic_rational_mul(&P->a, &I1->a, &I2->b);
        dyadic_rational_mul(&P->b, &I1->a, &I2->a);
      }
    }
  } else if (I2->is_point) {
    dyadic_interval_mul(P, I2, I1);
  } else {
    if (P->is_point) {
      dyadic_rational_construct(&P->b);
      P->is_point = 0;
    }

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //         = { x*y |

    dyadic_interval_t result;
    dyadic_interval_construct_zero(&result);

    dyadic_rational_t tmp;
    dyadic_rational_construct(&tmp);

    // I1.a x I2.a
    dyadic_rational_mul(&result.a, &I1->a, &I2->a);
    dyadic_rational_construct_copy(&result.b, &result.a); // was not constructed (is 0)
    result.a_open = result.b_open = I1->a_open || I2->a_open;
    result.is_point = 0;

    // I1.a x I2.b
    int tmp_open = I1->a_open || I2->b_open;
    dyadic_rational_mul(&tmp, &I1->a, &I2->b);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    tmp_open = I1->b_open || I2->a_open;
    dyadic_rational_mul(&tmp, &I1->b, &I2->a);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    tmp_open = I1->b_open || I2->b_open;
    dyadic_rational_mul(&tmp, &I1->b, &I2->b);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    dyadic_interval_swap(&result, P);
    dyadic_interval_destruct(&result);
    dyadic_rational_destruct(&tmp);
  }
}

int dyadic_interval_print(const dyadic_interval_t* I, FILE* out);

