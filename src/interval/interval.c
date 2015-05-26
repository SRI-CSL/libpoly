/*
 * interval.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include <rational_interval.h>
#include <dyadic_interval.h>
#include <interval.h>
#include <value.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

#include <assert.h>
#include <limits.h>

void lp_rational_interval_construct(lp_rational_interval_t* I,
    const lp_rational_t* a, int a_open,
    const lp_rational_t* b, int b_open)
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

void lp_dyadic_interval_construct(lp_dyadic_interval_t* I,
    const lp_dyadic_rational_t* a, int a_open,
    const lp_dyadic_rational_t* b, int b_open)
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

void lp_interval_construct(lp_interval_t* I,
    const lp_value_t* a, int a_open,
    const lp_value_t* b, int b_open)
{
  int cmp = lp_value_cmp(a, b);
  assert(cmp <= 0);
  lp_value_construct_copy(&I->a, a);
  if (cmp != 0) {
    lp_value_construct_copy(&I->b, b);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_rational_interval_construct_point(lp_rational_interval_t* I, const lp_rational_t* a)
{
  rational_construct_copy(&I->a, a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_construct_point(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
  dyadic_rational_construct_copy(&I->a, q);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_interval_construct_point(lp_interval_t* I, const lp_value_t* q) {
  lp_value_construct_copy(&I->a, q);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_rational_interval_construct_zero(lp_rational_interval_t* I) {
  rational_construct(&I->a);
  I->a_open = I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_construct_zero(lp_dyadic_interval_t* I) {
  dyadic_rational_construct(&I->a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_interval_construct_zero(lp_interval_t* I) {
  lp_value_construct_zero(&I->a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_rational_interval_construct_copy(lp_rational_interval_t* I, const lp_rational_interval_t* from) {
  rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_dyadic_interval_construct_copy(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from) {
  dyadic_rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    dyadic_rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_interval_construct_copy(lp_interval_t* I, const lp_interval_t* from) {
  lp_value_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    lp_value_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_interval_construct_full(lp_interval_t* I) {
  lp_value_construct(&I->a, LP_VALUE_MINUS_INFINITY, 0);
  lp_value_construct(&I->b, LP_VALUE_PLUS_INFINITY, 0);
  I->a_open = 1;
  I->b_open = 1;
  I->is_point = 0;
}

void lp_rational_interval_construct_from_dyadic(lp_rational_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open) {
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

void lp_rational_interval_construct_from_dyadic_interval(lp_rational_interval_t* I, const lp_dyadic_interval_t* from) {
  rational_construct_from_dyadic(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_from_dyadic(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_dyadic_interval_construct_from_dyadic(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_copy(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_rational_interval_construct_from_int(lp_rational_interval_t* I,
    long a, int a_open,
    long b, int b_open)
{
  assert(a <= b);
  rational_construct_from_int(&I->a, a, 0);
  if (a != b) {
    rational_construct_from_int(&I->b, b, 0);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_dyadic_interval_construct_from_int(lp_dyadic_interval_t* I,
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
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_rational_interval_construct_from_integer(lp_rational_interval_t* I,
    const lp_integer_t* a, int a_open,
    const lp_integer_t* b, int b_open)
{
  int cmp = integer_cmp(lp_Z, a, b);
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

void lp_dyadic_interval_construct_from_integer(lp_dyadic_interval_t* I,
    const lp_integer_t* a, int a_open,
    const lp_integer_t* b, int b_open)
{
  int cmp = integer_cmp(lp_Z, a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_from_integer(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_rational_interval_destruct(lp_rational_interval_t* I) {
  rational_destruct(&I->a);
  if (!I->is_point) {
    rational_destruct(&I->b);
  }
}

void lp_dyadic_interval_destruct(lp_dyadic_interval_t* I) {
  dyadic_rational_destruct(&I->a);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
}

void lp_interval_destruct(lp_interval_t* I) {
  lp_value_destruct(&I->a);
  if (!I->is_point) {
    lp_value_destruct(&I->b);
  }
}


void lp_rational_interval_assign(lp_rational_interval_t* I, const lp_rational_interval_t* from) {
  if (I != from) {
    if (I->is_point) {
      if (from->is_point) {
        rational_assign(&I->a, &from->a);
      } else {
        rational_assign(&I->a, &from->a);
        rational_construct_copy(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    } else {
      if (from->is_point) {
        rational_assign(&I->a, &from->a);
        rational_destruct(&I->b);
        I->a_open = I->b_open = 0;
        I->is_point = 1;
      } else {
        // both intervals
        rational_assign(&I->a, &from->a);
        rational_assign(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    }
  }
}

void lp_dyadic_interval_assign(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from) {
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

void lp_interval_assign(lp_interval_t* I, const lp_interval_t* from) {
  if (I != from) {
    if (I->is_point) {
      if (from->is_point) {
        lp_value_assign(&I->a, &from->a);
      } else {
        lp_value_assign(&I->a, &from->a);
        lp_value_construct_copy(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    } else {
      if (from->is_point) {
        lp_value_assign(&I->a, &from->a);
        lp_value_destruct(&I->b);
        I->a_open = I->b_open = 0;
        I->is_point = 1;
      } else {
        // both intervals
        lp_value_assign(&I->a, &from->a);
        lp_value_assign(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    }
  }
}

void lp_rational_interval_swap(lp_rational_interval_t* I1, lp_rational_interval_t* I2) {
  lp_rational_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

void lp_dyadic_interval_swap(lp_dyadic_interval_t* I1, lp_dyadic_interval_t* I2) {
  lp_dyadic_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2) {
  lp_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

int lp_rational_interval_sgn(const lp_rational_interval_t* I) {
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

int lp_dyadic_interval_sgn(const lp_dyadic_interval_t* I) {
  int a_sgn = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = dyadic_rational_sgn(&I->b);

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

int lp_dyadic_interval_cmp_integer(const lp_dyadic_interval_t* I, const lp_integer_t* z) {

  if (I->is_point) {
    return dyadic_rational_cmp_integer(&I->a, z);
  }

  // I = [a, b]

  int cmp_lower = dyadic_rational_cmp_integer(&I->a, z);
  if (cmp_lower > 0) {
    // a > z => [a, b] > z
    return 1;
  }
  if (cmp_lower == 0) {
    if (I->a_open) {
      // a == z => (a, b] > z
      return 1;
    } else {
      // a == z => [a, b] == z
      return 0;
    }
  }

  int cmp_upper = dyadic_rational_cmp_integer(&I->b, z);
  if (cmp_upper < 0) {
    // [a, b] < z
    return -1;
  }
  if (cmp_upper == 0) {
    if (I->b_open) {
      // [a, b) < z
      return -1;
    } else {
      // [a, b] == z
      return 0;
    }
  }

  // It's inside, return 0
  return 0;
}

int lp_dyadic_interval_cmp_dyadic_rational(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {

  if (I->is_point) {
    return dyadic_rational_cmp(&I->a, q);
  }

  // I = [a, b]

  int cmp_lower = dyadic_rational_cmp(&I->a, q);
  if (cmp_lower > 0) {
    // a > z => [a, b] > z
    return 1;
  }
  if (cmp_lower == 0) {
    if (I->a_open) {
      // a == z => (a, b] > z
      return 1;
    } else {
      // a == z => [a, b] == z
      return 0;
    }
  }

  int cmp_upper = dyadic_rational_cmp(&I->b, q);
  if (cmp_upper < 0) {
    // [a, b] < z
    return -1;
  }
  if (cmp_upper == 0) {
    if (I->b_open) {
      // [a, b) < z
      return -1;
    } else {
      // [a, b] == z
      return 0;
    }
  }

  // It's inside, return 0
  return 0;
}

int lp_dyadic_interval_cmp_rational(const lp_dyadic_interval_t* I, const lp_rational_t* q) {

  if (I->is_point) {
    return -rational_cmp_dyadic_rational(q, &I->a);
  }

  // I = [a, b]

  int cmp_a_q = -rational_cmp_dyadic_rational(q, &I->a);
  if (cmp_a_q > 0) {
    // a > q => [a, b] > q
    return 1;
  }
  if (cmp_a_q == 0) {
    if (I->a_open) {
      // a == q => (a, b] > q
      return 1;
    } else {
      // a == q => [a, b] == q
      return 0;
    }
  }

  int cmp_b_q = -rational_cmp_dyadic_rational(q, &I->b);
  if (cmp_b_q < 0) {
    // b < q => [a, b] < q
    return -1;
  }
  if (cmp_b_q == 0) {
    if (I->b_open) {
      // b == q => [a, b) < q
      return -1;
    } else {
      // b == q => [a, b] == q
      return 0;
    }
  }

  // It's inside, return 0
  return 0;
}


int lp_rational_interval_contains_integer(const lp_rational_interval_t* I, const lp_integer_t* z) {
  assert(0);
  (void)I;
  (void)z;
  return 1;
}

int lp_rational_interval_contains_dyadic_rational(const lp_rational_interval_t* I, const lp_dyadic_rational_t* dy_q) {
  assert(0);
  (void)I;
  (void)dy_q;
  return 1;
}

int lp_rational_interval_contains_algebraic_number(const lp_rational_interval_t* I, const lp_algebraic_number_t* a) {
  assert(0);
  (void)I;
  (void)a;
  return 1;
}

int lp_rational_interval_contains_value(const lp_rational_interval_t* I, const lp_value_t* v) {
  int cmp_a_v = -lp_value_cmp_rational(v, &I->a);
  if (I->is_point) {
    return cmp_a_v == 0;
  }
  if (I->a_open && cmp_a_v >= 0) return 0;
  if (!I->a_open && cmp_a_v > 0) return 0;
  int cmp_v_b = lp_value_cmp_rational(v, &I->b);
  if (I->b_open && cmp_v_b >= 0) return 0;
  if (!I->b_open && cmp_v_b > 0) return 0;
  return 1;
}

int lp_rational_interval_contains_rational(const lp_rational_interval_t* I, const lp_rational_t* q) {
  int cmp_a_q = rational_cmp(&I->a, q);
  if (I->is_point) {
    return cmp_a_q == 0;
  }
  if (I->a_open && cmp_a_q >= 0) return 0;
  if (!I->a_open && cmp_a_q > 0) return 0;
  int cmp_q_b = rational_cmp(q, &I->b);
  if (I->b_open && cmp_q_b >= 0) return 0;
  if (!I->b_open && cmp_q_b > 0) return 0;
  return 1;
}

int lp_interval_contains(const lp_interval_t* I, const lp_value_t* v) {
  int cmp_a_v = lp_value_cmp(&I->a, v);
  if (I->is_point) {
    return cmp_a_v == 0;
  }
  if (I->a_open && cmp_a_v >= 0) return 0;
  if (!I->a_open && cmp_a_v > 0) return 0;
  int cmp_v_b = lp_value_cmp(v, &I->b);
  if (I->b_open && cmp_v_b >= 0) return 0;
  if (!I->b_open && cmp_v_b > 0) return 0;
  return 1;
}


int lp_dyadic_interval_contains_dyadic_rational(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
  int cmp_a_q = dyadic_rational_cmp(&I->a, q);
  if (I->is_point) {
    return cmp_a_q == 0;
  }
  if (I->a_open && cmp_a_q >= 0) return 0;
  if (!I->a_open && cmp_a_q > 0) return 0;
  int cmp_q_b = dyadic_rational_cmp(q, &I->b);
  if (I->b_open && cmp_q_b >= 0) return 0;
  if (!I->b_open && cmp_q_b > 0) return 0;
  return 1;
}

int lp_rational_interval_contains_zero(const lp_rational_interval_t* I) {
  int sgn_a = rational_sgn(&I->a);
  if (I->is_point) {
    return sgn_a == 0;
  }
  if (I->a_open && sgn_a >= 0) return 0;
  if (!I->a_open && sgn_a > 0) return 0;
  int sgn_b = rational_sgn(&I->b);
  if (I->b_open && sgn_b <= 0) return 0;
  if (!I->b_open && sgn_b < 0) return 0;
  return 1;
}

int lp_dyadic_interval_contains_zero(const lp_dyadic_interval_t* I) {
  int sgn_a = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return sgn_a == 0;
  }
  if (I->a_open && sgn_a >= 0) return 0;
  if (!I->a_open && sgn_a > 0) return 0;
  int sgn_b = dyadic_rational_sgn(&I->b);
  if (I->b_open && sgn_b < 0) return 0;
  if (!I->b_open && sgn_b <= 0) return 0;
  return 1;
}


void lp_dyadic_interval_construct_from_split(lp_dyadic_interval_t* I_left, lp_dyadic_interval_t* I_right, const lp_dyadic_interval_t* I, int left_open, int right_open) {
  assert(!I->is_point);
  lp_dyadic_rational_t m;
  dyadic_rational_construct(&m);
  dyadic_rational_add(&m, &I->a, &I->b);
  dyadic_rational_div_2exp(&m, &m, 1);
  lp_dyadic_interval_construct(I_left, &I->a, I->a_open, &m, left_open);
  lp_dyadic_interval_construct(I_right, &m, right_open, &I->b, I->b_open);
  dyadic_rational_destruct(&m);
}

void lp_dyadic_interval_construct_intersection(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  if (I1->is_point) {
    assert(lp_dyadic_interval_contains_dyadic_rational(I2, &I1->a));
    lp_dyadic_interval_construct_copy(I, I1);
  } else if (I2->is_point) {
    assert(lp_dyadic_interval_contains_dyadic_rational(I1, &I2->a));
    lp_dyadic_interval_construct_copy(I, I2);
  } else {
    // (   [  )   ]
    int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
    const lp_dyadic_rational_t* max_a = cmp_a < 0 ? &I2->a : &I1->a;
    int a_open;
    if (cmp_a == 0) {
      a_open = I1->a_open || I2->a_open;
    } else if (cmp_a < 0) {
      a_open = I2->a_open;
    } else {
      a_open = I1->a_open;
    }

    int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
    const lp_dyadic_rational_t* min_b = cmp_b < 0 ? &I1->b : &I2->b;
    int b_open;
    if (cmp_b == 0) {
      b_open = I1->b_open || I2->b_open;
    } else if (cmp_b < 0) {
      b_open = I1->b_open;
    } else {
      b_open = I2->b_open;
    }

    lp_dyadic_interval_construct(I, max_a, a_open, min_b, b_open);
  }
}

void lp_dyadic_interval_collapse_to(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
  dyadic_rational_assign(&I->a, q);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_set_a(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(a, &I->b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->a, a);
    I->a_open = a_open;
  } else {
    assert(!a_open && !I->b_open);
    lp_dyadic_interval_collapse_to(I, a);
  }
}

void lp_dyadic_interval_set_b(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* b, int b_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(&I->a, b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->b, b);
    I->b_open = b_open;
  } else {
    assert(!I->a_open && !b_open);
    lp_dyadic_interval_collapse_to(I, b);
  }
}

int lp_dyadic_interval_equals(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
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

int lp_dyadic_interval_disjunct(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  if (I1->is_point) {
    return !lp_dyadic_interval_contains_dyadic_rational(I2, &I1->a);
  }
  if (I2->is_point) {
    return !lp_dyadic_interval_contains_dyadic_rational(I1, &I2->a);
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

void lp_dyadic_interval_scale(lp_dyadic_interval_t* I, int n) {
  assert(!I->is_point);
  if (n > 0) {
    dyadic_rational_mul_2exp(&I->a, &I->a, n);
    dyadic_rational_mul_2exp(&I->b, &I->b, n);
  } else {
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
  }
}

/** Prints the interval to the given stream. */
int lp_rational_interval_print(const lp_rational_interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += rational_print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += rational_print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += rational_print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

int lp_dyadic_interval_print(const lp_dyadic_interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += dyadic_rational_print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += dyadic_rational_print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += dyadic_rational_print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

int lp_interval_print(const lp_interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += lp_value_print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += lp_value_print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += lp_value_print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

int lp_dyadic_interval_is_point(const lp_dyadic_interval_t* I) {
  return I->is_point;
}

int lp_rational_interval_is_point(const lp_rational_interval_t* I) {
  return I->is_point;
}

int lp_interval_is_point(const lp_interval_t* I) {
  return I->is_point;
}

const lp_dyadic_rational_t* lp_dyadic_interval_get_point(const lp_dyadic_interval_t* I) {
  return &I->a;
}

const lp_rational_t* lp_rational_interval_get_point(const lp_rational_interval_t* I) {
  return &I->a;
}

const lp_value_t* lp_interval_get_point(const lp_interval_t* I) {
  return &I->a;
}

const lp_value_t* lp_interval_get_lower_bound(const lp_interval_t* I) {
  return &I->a;
}

const lp_value_t* lp_interval_get_upper_bound(const lp_interval_t* I) {
  return I->is_point ? &I->a : &I->b;
}

void lp_interval_pick_value(const lp_interval_t* I, lp_value_t* v) {
  if (I->is_point) {
    lp_value_assign(v, &I->a);
  } else {
    lp_value_get_value_between(&I->a, I->a_open, &I->b, I->b_open, v);
  }
}

int lp_dyadic_interval_size(const lp_dyadic_interval_t* I) {
  // If point, size is 0
  if (I->is_point) {
    return 0;
  } else {
    return dyadic_rational_get_distance_size(&I->a, &I->b);
  }
}

int lp_interval_size_approx(const lp_interval_t* I) {
  if (I->is_point) {
    return INT_MIN;
  } else {
    return lp_value_get_distance_size_approx(&I->a, &I->b);
  }
}

char* lp_interval_to_string(const lp_interval_t* I) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_interval_print(I, f);
  fclose(f);
  return str;
}

int lp_interval_cmp_lower_bounds(const lp_interval_t* I1, const lp_interval_t* I2) {
  const lp_value_t* I1_lb = lp_interval_get_lower_bound(I1);
  const lp_value_t* I2_lb = lp_interval_get_lower_bound(I2);

  int cmp = lp_value_cmp(I1_lb, I2_lb);
  if (cmp != 0) {
    return cmp;
  } else {
    // Equal, so we check if they are strict
    if (I1->a_open == I2->a_open) {
      return 0;
    } else {
      if (I1->a_open) {
        // (a, b) > [a, b)
        return 1;
      } else {
        // [a, b) < (a, b)
        return -1;
      }
    }
  }
}

int lp_interval_cmp_upper_bounds(const lp_interval_t* I1, const lp_interval_t* I2) {
  const lp_value_t* I1_ub = lp_interval_get_upper_bound(I1);
  const lp_value_t* I2_ub = lp_interval_get_upper_bound(I2);

  int cmp = lp_value_cmp(I1_ub, I2_ub);
  if (cmp != 0) {
    return cmp;
  } else {
    // Equal, so we check if they are strict
    if (I1->b_open == I2->b_open) {
      return 0;
    } else {
      if (I1->b_open) {
        // (a, b) < (a, b]
        return -1;
      } else {
        // (a, b] > (a, b)
        return 1;
      }
    }
  }
}

lp_interval_cmp_t lp_interval_cmp(const lp_interval_t* I1, const lp_interval_t* I2) {
  return lp_interval_cmp_with_intersect(I1, I2, 0);
}

lp_interval_cmp_t lp_interval_cmp_with_intersect(const lp_interval_t* I1, const lp_interval_t* I2, lp_interval_t* P) {

  int cmp_ub = lp_interval_cmp_upper_bounds(I1, I2);
  int cmp_lb = lp_interval_cmp_lower_bounds(I1, I2);

  if (cmp_ub == 0 && cmp_lb == 0) {
    if (P) {
      lp_interval_t result;
      lp_interval_construct_copy(&result, I1);
      lp_interval_swap(&result, P);
      lp_interval_destruct(&result);
    }
    return LP_INTERVAL_CMP_EQ;
  }

  if (cmp_ub < 0 && cmp_lb > 0) {
    //   ( )
    // (      )
    if (P) {
      lp_interval_t result;
      lp_interval_construct_copy(&result, I1);
      lp_interval_swap(&result, P);
      lp_interval_destruct(&result);
    }
    return LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1;
  }

  if (cmp_ub > 0 && cmp_lb < 0) {
    // (       )
    //    ( )
    if (P) {
      lp_interval_t result;
      lp_interval_construct_copy(&result, I2);
      lp_interval_swap(&result, P);
      lp_interval_destruct(&result);
    }
    return LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2;
  }

  if (cmp_ub == 0) {
    if (cmp_lb > 0) {
      //   (    )
      // (      )
      if (P) {
        lp_interval_t result;
        lp_interval_construct_copy(&result, I1);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_GEQ_WITH_INTERSECT_I1;
    }
    if (cmp_lb < 0) {
      // (    )
      //    ( )
      if (P) {
        lp_interval_t result;
        lp_interval_construct_copy(&result, I2);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_LEQ_WITH_INTERSECT_I2;
    }
  }

  if (cmp_lb == 0) {
    if (cmp_ub > 0) {
      // (     )
      // (   )
      if (P) {
        lp_interval_t result;
        lp_interval_construct_copy(&result, I2);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2;
    }
    if (cmp_ub < 0) {
      // (  )
      // (     )
      if (P) {
        lp_interval_t result;
        lp_interval_construct_copy(&result, I1);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1;
    }
  }

  /**
   * Here we know that both comparisons go the same way
   *
   *  (   )        (    )
   *    (   )    (    )
   *
   *  ( )            ( )
   *      ( )    ( )
   */

  if (cmp_ub < 0) {
    assert(cmp_lb < 0);
    const lp_value_t* I1_ub = lp_interval_get_upper_bound(I1);
    const lp_value_t* I2_lb = lp_interval_get_lower_bound(I2);
    int cmp_I1_ub_I2_lb = lp_value_cmp(I1_ub, I2_lb);
    if (cmp_I1_ub_I2_lb == 0 && (I1->b_open || I2->a_open)) {
      cmp_I1_ub_I2_lb = -1;
    }
    if (cmp_I1_ub_I2_lb == 0) {
      // I1: (  ]
      // I2:    [  )
      assert(!I1->b_open && !I2->a_open);
      if (P) {
        lp_interval_t result;
        lp_interval_construct_point(&result, &I2->a);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_LT_WITH_INTERSECT;
    } else if (cmp_I1_ub_I2_lb < 0) {
      // I1: (  )
      // I2      (  )
      return LP_INTERVAL_CMP_LT_NO_INTERSECT;
    } else {
      // I1: (   )
      // I2:   (   )
      if (P) {
        lp_interval_t result;
        lp_interval_construct(&result, I2_lb, I2->a_open, I1_ub, I1->b_open);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_LT_WITH_INTERSECT;
    }
  } else {
    assert(cmp_ub > 0);
    assert(cmp_lb > 0);
    const lp_value_t* I1_lb = lp_interval_get_lower_bound(I1);
    const lp_value_t* I2_ub = lp_interval_get_upper_bound(I2);
    int cmp_I1_lb_I2_ub = lp_value_cmp(I1_lb, I2_ub);
    if (cmp_I1_lb_I2_ub == 0 && (I1->a_open || I2->b_open)) {
      cmp_I1_lb_I2_ub = -1;
    }
    if (cmp_I1_lb_I2_ub == 0) {
      // I1:    [  )
      // I2: (  ]
      assert(!I1->a_open && !I2->b_open);
      if (P) {
        lp_interval_t result;
        lp_interval_construct_point(&result, &I1->a);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_GT_WITH_INTERSECT;
    } else if (cmp_I1_lb_I2_ub < 0) {
      // I1:   (   )
      // I2: (   )
      if (P) {
        lp_interval_t result;
        lp_interval_construct(&result, I1_lb, I1->a_open, I2_ub, I2->b_open);
        lp_interval_swap(&result, P);
        lp_interval_destruct(&result);
      }
      return LP_INTERVAL_CMP_GT_WITH_INTERSECT;
    } else {
      // I1:     (  )
      // I2: (  )
      return LP_INTERVAL_CMP_GT_NO_INTERSECT;
    }
  }
}
