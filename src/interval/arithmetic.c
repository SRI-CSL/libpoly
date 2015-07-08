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

#include "interval/arithmetic.h"

#include "number/rational.h"
#include "number/dyadic_rational.h"

int rational_interval_endpoint_lt(const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open) {
  int cmp = rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

int dyadic_interval_endpoint_lt(const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t*b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

void rational_interval_add(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {

  if (I1->is_point && I2->is_point) {
    if (!S->is_point) {
      rational_destruct(&S->b);
    }
    rational_add(&S->a, &I1->a, &I2->a);
    S->b_open = S->a_open = 0;
    S->is_point = 1;
    return;
  }

  if (I2->is_point) {
    // Reuse symmetry
    rational_interval_add(S, I2, I1);
    return;
  }

  lp_rational_interval_t result;

  if (I1->is_point) {
    // Just shift by I1->a (I2 is not a point)
    lp_rational_interval_construct_copy(&result, I2);
    rational_add(&result.a, &result.a, &I1->a);
    rational_add(&result.b, &result.b, &I1->a);
  } else {
    // [a, b] + [c, d] = [a + c, b + d]
    rational_construct(&result.a);
    rational_construct(&result.b);
    rational_add(&result.a, &I1->a, &I2->a);
    rational_add(&result.b, &I1->b, &I2->b);
    result.a_open = I1->a_open || I2->a_open;
    result.b_open = I1->b_open || I2->b_open;
    result.is_point = 0;
  }

  lp_rational_interval_swap(&result, S);
  lp_rational_interval_destruct(&result);
}

void dyadic_interval_add(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {

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

  lp_dyadic_interval_t result;

  if (I1->is_point) {
    // Just shift by I1->a
    lp_dyadic_interval_construct_copy(&result, I2);
    dyadic_rational_add(&result.a, &result.a, &I1->a);
    dyadic_rational_add(&result.b, &result.b, &I1->a);
    return;
  } else {
    // Both non-points
    // [a, b] + [c, d] = [a + c, b + d]
    dyadic_rational_construct(&result.a);
    dyadic_rational_construct(&result.b);
    dyadic_rational_add(&result.a, &I1->a, &I2->a);
    dyadic_rational_add(&result.b, &I1->b, &I2->b);
    result.a_open = I1->a_open || I2->a_open;
    result.b_open = I1->b_open || I2->b_open;
    result.is_point = 0;
  }

  lp_dyadic_interval_swap(&result, S);
  lp_dyadic_interval_destruct(&result);
}

void rational_interval_neg(lp_rational_interval_t* N, const lp_rational_interval_t* I) {
  if (I->is_point) {
    if (!N->is_point) {
      rational_destruct(&N->b);
    }
    rational_neg(&N->a, &I->a);
    N->b_open = N->a_open = 0;
    N->is_point = 1;
    return;
  }

  if (N->is_point) {
    rational_construct(&N->b);
    N->is_point = 0;
  }

  // -[a, b] = [-b, -a]
  // (doing swap in case I and N are the same)
  rational_neg(&N->a, &I->a);
  rational_neg(&N->b, &I->b);
  N->a_open = I->a_open;
  N->b_open = I->b_open;

  rational_swap(&N->a, &N->b);
  size_t tmp = N->a_open;
  N->a_open = N->b_open;
  N->b_open = tmp;
}

void dyadic_interval_neg(lp_dyadic_interval_t* N, const lp_dyadic_interval_t* I) {
  if (I->is_point) {
    if (!N->is_point) {
      dyadic_rational_destruct(&N->b);
    }
    dyadic_rational_neg(&N->a, &I->a);
    N->b_open = N->a_open = 0;
    N->is_point = 1;
    return;
  }

  if (N->is_point) {
    dyadic_rational_construct(&N->b);
    N->is_point = 0;
  }

  // -[a, b] = [-b, -a]
  // (doing swap in case I and N are the same)
  dyadic_rational_neg(&N->a, &I->a);
  dyadic_rational_neg(&N->b, &I->b);
  N->a_open = I->a_open;
  N->b_open = I->b_open;

  dyadic_rational_swap(&N->a, &N->b);
  size_t tmp = N->a_open;
  N->a_open = N->b_open;
  N->b_open = tmp;
}

void rational_interval_sub(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {
  lp_rational_interval_t neg;
  lp_rational_interval_construct_copy(&neg, I2);
  rational_interval_neg(&neg, &neg);
  rational_interval_add(S, I1, &neg);
  lp_rational_interval_destruct(&neg);
}

void dyadic_interval_sub(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  lp_dyadic_interval_t neg;
  lp_dyadic_interval_construct_copy(&neg, I2);
  dyadic_interval_neg(&neg, &neg);
  dyadic_interval_add(S, I1, &neg);
  lp_dyadic_interval_destruct(&neg);
}

void rational_interval_mul(lp_rational_interval_t* P, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {
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
      // I1 is point
      // I2 is not point

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
        lp_rational_interval_t result;
        rational_construct(&result.a);
        rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        rational_mul(&result.a, &I1->a, &I2->a);
        rational_mul(&result.b, &I1->a, &I2->b);
        lp_rational_interval_swap(&result, P);
        lp_rational_interval_destruct(&result);
      } else {
        lp_rational_interval_t result;
        rational_construct(&result.a);
        rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        rational_mul(&result.b, &I1->a, &I2->a);
        rational_mul(&result.a, &I1->a, &I2->b);
        lp_rational_interval_swap(&result, P);
        lp_rational_interval_destruct(&result);
      }
    }
  } else if (I2->is_point) {
    rational_interval_mul(P, I2, I1);
  } else {
    if (P->is_point) {
      rational_construct(&P->b);
      P->is_point = 0;
    }

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //         = { x*y |

    lp_rational_interval_t result;
    lp_rational_interval_construct_zero(&result);

    lp_rational_t tmp;
    rational_construct(&tmp);

    // I1.a x I2.a
    rational_mul(&result.a, &I1->a, &I2->a);
    rational_construct_copy(&result.b, &result.a);
    result.a_open = result.b_open = I1->a_open || I2->a_open;
    result.is_point = 0;

    // I1.a x I2.b
    int tmp_open = I1->a_open || I2->b_open;
    rational_mul(&tmp, &I1->a, &I2->b);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    tmp_open = I1->b_open || I2->a_open;
    rational_mul(&tmp, &I1->b, &I2->a);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    tmp_open = I1->b_open || I2->b_open;
    rational_mul(&tmp, &I1->b, &I2->b);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    lp_rational_interval_swap(&result, P);
    lp_rational_interval_destruct(&result);
    rational_destruct(&tmp);
  }
}

void dyadic_interval_mul(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
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
        lp_dyadic_interval_t result;
        dyadic_rational_construct(&result.a);
        dyadic_rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        dyadic_rational_mul(&result.a, &I1->a, &I2->a);
        dyadic_rational_mul(&result.b, &I1->a, &I2->b);
        lp_dyadic_interval_swap(&result, P);
        lp_dyadic_interval_destruct(&result);
      } else {
        lp_dyadic_interval_t result;
        dyadic_rational_construct(&result.a);
        dyadic_rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->b_open;
        result.b_open = I2->a_open;
        dyadic_rational_mul(&result.a, &I1->a, &I2->b);
        dyadic_rational_mul(&result.b, &I1->a, &I2->a);
        lp_dyadic_interval_swap(&result, P);
        lp_dyadic_interval_destruct(&result);
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

    lp_dyadic_interval_t result;
    lp_dyadic_interval_construct_zero(&result);

    lp_dyadic_rational_t tmp;
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

    lp_dyadic_interval_swap(&result, P);
    lp_dyadic_interval_destruct(&result);
    dyadic_rational_destruct(&tmp);
  }
}

void rational_interval_pow(lp_rational_interval_t* P, const lp_rational_interval_t* I, unsigned n) {
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
      int sgn = lp_rational_interval_sgn(I);
      rational_pow(&P->a, &I->a, n);
      rational_pow(&P->b, &I->b, n);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        if (rational_interval_endpoint_lt(&P->b, I->b_open, &P->a, I->a_open)) {
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

void dyadic_interval_pow(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I, unsigned n) {
  if (n == 0) {
    // I^0 = [1]
    if (!P->is_point) {
      P->is_point = 1;
      dyadic_rational_destruct(&P->b);
    }
    dyadic_rational_assign_int(&P->a, 1, 1);
    P->a_open = 0;
    P->b_open = 0;
  } else if (I->is_point) {
    // Plain power
    if (!P->is_point) {
      dyadic_rational_destruct(&P->b);
      P->is_point = 1;
      P->a_open = P->b_open = 0;
    }
    dyadic_rational_pow(&P->a, &I->a, n);
  } else {
    if (P->is_point) {
      P->is_point = 0;
      dyadic_rational_construct(&P->b);
    }
    if (n % 2) {
      // For odd powers we are monotonic, i.e. [a, b]^n = [a^n, b^n]
      P->a_open = I->a_open;
      P->b_open = I->b_open;
      dyadic_rational_pow(&P->a, &I->a, n);
      dyadic_rational_pow(&P->b, &I->b, n);
    } else {
      // Even powers depend on whether 0 is in the interval
      int sgn = lp_dyadic_interval_sgn(I);
      dyadic_rational_pow(&P->a, &I->a, n);
      dyadic_rational_pow(&P->b, &I->b, n);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        if (dyadic_interval_endpoint_lt(&P->b, I->b_open, &P->a, I->a_open)) {
          dyadic_rational_swap(&P->b, &P->a);
          P->b_open = I->a_open;
        }
        dyadic_rational_assign_int(&P->a, 0, 1);
        P->a_open = 0;
      } else if (sgn > 0) {
        // P = I^n
        P->a_open = I->a_open;
        P->b_open = I->b_open;
      } else {
        // negative turns positive, so we flip
        dyadic_rational_swap(&P->a, &P->b);
        P->a_open = I->b_open;
        P->b_open = I->a_open;
      }
    }
  }
}
