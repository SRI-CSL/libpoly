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
#include <monomial.h>

#include "polynomial/gcd.h"
#include "polynomial/output.h"
#include "polynomial/polynomial_vector.h"

#include "upolynomial/upolynomial.h"

#include "utils/statistics.h"
#include "utils/debug_trace.h"

void monomial_gcd_visit(const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* data) {
  lp_monomial_t* gcd = (lp_monomial_t*) data;
  if (integer_is_zero(ctx->K, &gcd->a)) {
    lp_monomial_assign(ctx, gcd, m, 0);
  } else {
    lp_monomial_gcd(ctx, gcd, gcd, m);
  }
}

/**
 * Extracts the largest monomial power out of P and Q and into gcd, also divide.
 * For example, P and Q in Z[y, x]
 *
 *  P = 4*y*x^2 + 2*y^2 = 2*y^2*(2*x^2 + 1)
 *  Q = 2*y^3*x^3
 *
 * gives
 *
 *  gcd = 2*y^2
 *  P = 2*x^2 + 1
 *  Q = 2*y*x^3
 */
void coefficient_gcd_monomial_extract(const lp_polynomial_context_t* ctx, coefficient_t* gcd, coefficient_t* P, coefficient_t* Q) {

  TRACE("coefficient", "coefficient_gcd_monomial_extract()\n");

  if (trace_is_enabled("coefficient")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
  }

  assert(P != Q);

  lp_monomial_t m_P_gcd, m_Q_gcd, m_tmp;
  lp_monomial_construct(ctx, &m_P_gcd);
  lp_monomial_construct(ctx, &m_Q_gcd);
  lp_monomial_construct(ctx, &m_tmp);

  // Compute the gcd
  coefficient_traverse(ctx, P, monomial_gcd_visit, &m_tmp, &m_P_gcd);
  lp_monomial_clear(ctx, &m_tmp);
  coefficient_traverse(ctx, Q, monomial_gcd_visit, &m_tmp, &m_Q_gcd);

  if (trace_is_enabled("coefficient")) {
    tracef("P_gcd = "); monomial_print(ctx, &m_P_gcd, trace_out); tracef("\n");
    tracef("Q_gcd = "); monomial_print(ctx, &m_Q_gcd, trace_out); tracef("\n");
  }

  // Final gcd
  lp_monomial_t m_gcd;
  lp_monomial_construct(ctx, &m_gcd);
  lp_monomial_gcd(ctx, &m_gcd, &m_P_gcd, &m_Q_gcd);

  // Construct the result
  coefficient_t result;
  coefficient_construct(ctx, &result);
  coefficient_add_ordered_monomial(ctx, &m_gcd, &result);

  // Divide P and Q with their gcds
  coefficient_t P_gcd, Q_gcd;
  coefficient_construct(ctx, &P_gcd);
  coefficient_construct(ctx, &Q_gcd);
  coefficient_add_ordered_monomial(ctx, &m_P_gcd, &P_gcd);
  coefficient_add_ordered_monomial(ctx, &m_Q_gcd, &Q_gcd);
  coefficient_div(ctx, P, P, &P_gcd);
  coefficient_div(ctx, Q, Q, &Q_gcd);
  coefficient_destruct(&P_gcd);
  coefficient_destruct(&Q_gcd);

  // Output the result
  coefficient_swap(&result, gcd);
  coefficient_destruct(&result);

  lp_monomial_destruct(&m_gcd);
  lp_monomial_destruct(&m_tmp);
  lp_monomial_destruct(&m_Q_gcd);
  lp_monomial_destruct(&m_P_gcd);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd_monomial_extract() =>"); coefficient_print(ctx, gcd, trace_out); tracef("\n");
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
  }
}

/**
 * Takes two (primitive) coefficients over the same variable, makes them univariate by
 * substituting 0 for other variables (if any). Then it computes the
 * univariate gcd of these. If the coefficients were univariate already, or
 * the result is a constant (i.e. gcd = 1), the result is precise.
 */
int coefficient_gcd_pp_univariate(const lp_polynomial_context_t* ctx,
    coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2) {

  assert(C1->type == COEFFICIENT_POLYNOMIAL);
  assert(C2->type == COEFFICIENT_POLYNOMIAL);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp_univariate()\n");
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
  }

  int C1_vanishes = integer_is_zero(ctx->K, coefficient_get_constant(coefficient_lc(C1)));
  int C2_vanishes = integer_is_zero(ctx->K, coefficient_get_constant(coefficient_lc(C2)));

  if (C1_vanishes || C2_vanishes) {
    // One of C1 or C2 vanishes in the univariate conversion, we're not precise enough
    return 0;
  }

  lp_variable_t x = VAR(C1);
  assert(x == VAR(C2));

  lp_upolynomial_t* C1_u = coefficient_to_univariate(ctx, C1);
  lp_upolynomial_t* C2_u = coefficient_to_univariate(ctx, C2);
  lp_upolynomial_t* gcd_u = lp_upolynomial_gcd(C1_u, C2_u);

  coefficient_t gcd_tmp;
  coefficient_construct_from_univariate(ctx, &gcd_tmp, gcd_u, x);
  coefficient_swap(&gcd_tmp, gcd);
  coefficient_destruct(&gcd_tmp);

  lp_upolynomial_delete(C1_u);
  lp_upolynomial_delete(C2_u);
  lp_upolynomial_delete(gcd_u);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp_univariate() => ");
    tracef("gcd = "); coefficient_print(ctx, gcd, trace_out); tracef("\n");
  }

  if (gcd->type == COEFFICIENT_NUMERIC) {
    integer_assign_int(ctx->K, &gcd->value.num, 1);
    return 1;
  } else if (coefficient_is_univariate(C1) && coefficient_is_univariate(C2)) {
    return 1;
  } else {
    return 0;
  }
}

STAT_DECLARE(int, coefficient, gcd_pp_euclid)

/**
 * Compute the gcd of two primitive polynomials P and Q. The polynomials P and
 * Q will be used and changed in the computation.
 */
void coefficient_gcd_pp_euclid(const lp_polynomial_context_t* ctx, coefficient_t* gcd, coefficient_t* P, coefficient_t* Q) {

  TRACE("coefficient", "coefficient_gcd_pp()\n");
  STAT(coefficient, gcd_pp_euclid) ++;

  if (trace_is_enabled("coefficient::gcd")) {
    tracef("gcd\n")
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
  }

  // Try to compute the univariate GCD first
  coefficient_t gcd_u;
  coefficient_construct(ctx, &gcd_u);

  int precise = coefficient_gcd_pp_univariate(ctx, &gcd_u, P, Q);
  if (precise) {
    // GCD = 1, just copy the univariate gcd
    coefficient_swap(gcd, &gcd_u);
  } else {

    coefficient_t R;
    coefficient_construct(ctx, &R);

    //
    // We compute the reduction of P and Q in Z[y, x], i.e.
    //
    //   a*P = b*Q + R
    //
    // with a in Z[y], b in Z[y, x], and deg(R) < deg(Q) or deg(R) == 0.
    //
    // P and Q are primitive so GCD(P, Q) should be primitive, i.e. in
    // Z[y, x] or 1. therefore GCD(P, Q) = 1 if R != 0, or
    // GCD(P, Q) = po(Q) if R = 0
    //
    do {

      // One step reduction
      coefficient_reduce(ctx, P, Q, 0, 0, &R, REMAINDERING_PSEUDO_SPARSE);

      int cmp_type = coefficient_cmp_type(ctx, Q, &R);
      if (cmp_type == 0) {
        // P = Q
        // Q = pp(R)
        coefficient_swap(P, Q);
        coefficient_swap(Q, &R);
        coefficient_pp(ctx, Q, Q);
      } else {
        assert(cmp_type > 0);
        if (!coefficient_is_zero(ctx, &R)) {
          coefficient_destruct(Q);
          coefficient_construct_from_int(ctx, Q, 1);
        }
        break;
      }
    } while (1);

    coefficient_swap(Q, gcd);
    coefficient_destruct(&R);
  }

  coefficient_destruct(&gcd_u);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp() => "); coefficient_print(ctx, gcd, trace_out); tracef("\n");
  }
}

STAT_DECLARE(int, coefficient, gcd_pp_subresultant)

/**
 * Compute the gcd of two primitive polynomials P and Q. The polynomials P and
 * Q will be used and changed in the computation.
 */
void coefficient_gcd_pp_subresultant(const lp_polynomial_context_t* ctx, coefficient_t* gcd, coefficient_t* P, coefficient_t* Q) {

  TRACE("coefficient", "coefficient_gcd_pp_euclid()\n");
  STAT(coefficient, gcd_pp_subresultant) ++;

  if (trace_is_enabled("coefficient::gcd")) {
    tracef("gcd\n")
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
  }

  // Try to compute the univariate GCD first
  coefficient_t gcd_u;
  coefficient_construct(ctx, &gcd_u);

  int precise = coefficient_gcd_pp_univariate(ctx, &gcd_u, P, Q);
  if (precise) {
    // GCD = 1, just copy the univariate gcd
    coefficient_swap(gcd, &gcd_u);
  } else {

    // Make sure that P >= Q
    if (SIZE(P) < SIZE(Q)) {
      coefficient_t* tmp = P; P = Q; Q = tmp;
    }

    coefficient_t R;
    coefficient_construct(ctx, &R);

    coefficient_t h, g;
    coefficient_construct_from_int(ctx, &g, 1);
    coefficient_construct_from_int(ctx, &h, 1);

    coefficient_t tmp1, tmp2;
    coefficient_construct(ctx, &tmp1);
    coefficient_construct(ctx, &tmp2);

    // Subresultant GCD
    //
    do {

      // d = deg(P) - deg(Q)
      assert(SIZE(P) >= SIZE(Q));
      unsigned delta = SIZE(P) - SIZE(Q);

      // One step reduction
      coefficient_reduce(ctx, P, Q, 0, 0, &R, REMAINDERING_PSEUDO_SPARSE);

      if (trace_is_enabled("coefficient::gcd")) {
        tracef("------------\n");
        tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
        tracef("Q = "); coefficient_print(ctx, Q, trace_out); tracef("\n");
        tracef("h = "); coefficient_print(ctx, &h, trace_out); tracef("\n");
        tracef("g = "); coefficient_print(ctx, &g, trace_out); tracef("\n");
        tracef("d = %u\n", delta);
      }

      int cmp_type = coefficient_cmp_type(ctx, Q, &R);
      if (cmp_type == 0) {
        // P = Q
        coefficient_swap(P, Q);
        // Q = R/g*(h^delta)
        coefficient_div(ctx, &tmp1, &R, &g);
        coefficient_pow(ctx, &tmp2, &h, delta);
        coefficient_div(ctx, Q, &tmp1, &tmp2);
        // g = lc(P)
        coefficient_assign(ctx, &g, coefficient_lc(P));
        // h = h^(1-delta)*g^delta
        if (delta == 0) {
          // h = h, nothing to do
        } else if (delta == 1) {
          // h = g
          coefficient_assign(ctx, &h, &g);
        } else {
          // h = g^delta/h^(delta-1))
          coefficient_pow(ctx, &tmp1, &g, delta);
          coefficient_pow(ctx, &tmp2, &h, delta-1);
          coefficient_div(ctx, &h, &tmp1, &tmp2);
        }
      } else {
        assert(cmp_type > 0);
        if (!coefficient_is_zero(ctx, &R)) {
          coefficient_destruct(Q);
          coefficient_construct_from_int(ctx, Q, 1);
        } else {
          coefficient_pp(ctx, Q, Q);
        }
        break;
      }
    } while (1);

    coefficient_swap(Q, gcd);
    coefficient_destruct(&R);

    coefficient_destruct(&h);
    coefficient_destruct(&g);
    coefficient_destruct(&tmp1);
    coefficient_destruct(&tmp2);
  }

  coefficient_destruct(&gcd_u);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp() => "); coefficient_print(ctx, gcd, trace_out); tracef("\n");
  }
}

STAT_DECLARE(int, coefficient, gcd)

void coefficient_gcd(const lp_polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_gcd()\n");
  STAT(coefficient, gcd) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
  }

  assert(ctx->K == lp_Z);

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  if (cmp_type < 0) {
    const coefficient_t* tmp = C1;
    C1 = C2;
    C2 = tmp;
    cmp_type = -cmp_type;
  }

  if (cmp_type == 0) {
    switch (C1->type) {
    case COEFFICIENT_NUMERIC:
      if (gcd->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(gcd);
        coefficient_construct(ctx, gcd);
      }
      integer_gcd_Z(&gcd->value.num, &C1->value.num, &C2->value.num);
      break;
    case COEFFICIENT_POLYNOMIAL:
    {
      coefficient_t P, Q;
      if (SIZE(C1) > SIZE(C2)) {
        coefficient_construct_copy(ctx, &P, C1);
        coefficient_construct_copy(ctx, &Q, C2);
      } else {
        coefficient_construct_copy(ctx, &P, C2);
        coefficient_construct_copy(ctx, &Q, C1);
      }

      // Get the common power variables out
      coefficient_t gcd_mon;
      coefficient_construct(ctx, &gcd_mon);
      coefficient_gcd_monomial_extract(ctx, &gcd_mon, &P, &Q);

      // If monomial extraction changed the type, we need to go again
      if (coefficient_cmp_type(ctx, C1, &P) != 0 || coefficient_cmp_type(ctx, C2, &Q) != 0) {
        coefficient_gcd(ctx, gcd, &P, &Q);
      } else {
        // Normalize the P and Q to be primitive (and keep the content)
        coefficient_t P_cont, Q_cont;
        coefficient_construct(ctx, &P_cont);
        coefficient_construct(ctx, &Q_cont);
        coefficient_pp_cont(ctx, &P, &P_cont, &P);
        coefficient_pp_cont(ctx, &Q, &Q_cont, &Q);

        // Get the gcd of the content
        coefficient_t gcd_cont;
        coefficient_construct(ctx, &gcd_cont);
        coefficient_gcd(ctx, &gcd_cont, &P_cont, &Q_cont);

        // Get the gcd of the primitive parts
        coefficient_gcd_pp_euclid(ctx, gcd, &P, &Q);

        // Multiply in the content gcd
        coefficient_mul(ctx, gcd, gcd, &gcd_cont);

        coefficient_destruct(&P_cont);
        coefficient_destruct(&Q_cont);
        coefficient_destruct(&gcd_cont);
      }

      // Multiply in the monomial gcd
      coefficient_mul(ctx, gcd, gcd, &gcd_mon);

      // Remove temps
      coefficient_destruct(&P);
      coefficient_destruct(&Q);
      coefficient_destruct(&gcd_mon);
      break;
    }
    default:
      assert(0);
      break;
    }
  } else {
    // C1 in Z[y, x]
    // C2 in Z[y]
    // so GCD(C1, C2) = GCD(cont(C1), C2)
    coefficient_t cont;
    coefficient_construct(ctx, &cont);
    coefficient_cont(ctx, &cont, C1);
    coefficient_gcd(ctx, gcd, &cont, C2);
    coefficient_destruct(&cont);
  }

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_gcd() => "); coefficient_print(ctx, gcd, trace_out); tracef("\n");
  }

  if (trace_is_enabled("coefficient::gcd::sage")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
    tracef("gcd = "); coefficient_print(ctx, gcd, trace_out); tracef("\n");
    tracef("gcd_sage = C1.gcd(C2)\n");
    tracef("if (gcd != gcd_sage):\n");
    tracef("\tprint 'C1 =', C1\n");
    tracef("\tprint 'C2 =', C2\n");
  }

  assert(coefficient_is_normalized(ctx, gcd));
}

STAT_DECLARE(int, coefficient, lcm)

void coefficient_lcm(const lp_polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_lcm()\n");
  STAT(coefficient, lcm) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
  }

  assert(ctx->K == lp_Z);

  if (C1->type == COEFFICIENT_NUMERIC && C2->type == COEFFICIENT_NUMERIC) {
    // Integer LCM
    if (lcm->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(lcm);
      coefficient_construct(ctx, lcm);
    }
    integer_lcm_Z(&lcm->value.num, &C1->value.num, &C2->value.num);
  } else {
    // LCM(C1, C2) = C1*C2/GCD(C1, C2)
    coefficient_t gcd;
    coefficient_construct(ctx, &gcd);
    coefficient_gcd(ctx, &gcd, C1, C2);
    if (coefficient_is_one(ctx, &gcd)) {
      coefficient_mul(ctx, lcm, C1, C2);
    } else {
      if (coefficient_cmp_type(ctx, C1, C2) <= 0) {
        coefficient_div(ctx, lcm, C1, &gcd);
        coefficient_mul(ctx, lcm, lcm, C2);
      } else {
        coefficient_div(ctx, lcm, C2, &gcd);
        coefficient_mul(ctx, lcm, lcm, C1);
      }
    }
    if (coefficient_lc_sgn(ctx, lcm) < 0) {
      coefficient_neg(ctx, lcm, lcm);
    }
    coefficient_destruct(&gcd);
  }

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_lcm() => "); coefficient_print(ctx, lcm, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, lcm));
}

STAT_DECLARE(int, coefficient, pp_cont)

void coefficient_pp_cont(const lp_polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C) {

  TRACE("coefficient", "coefficient_pp_cont()\n");
  STAT(coefficient, pp_cont) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
  }

  assert(ctx->K == lp_Z);

  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    if (cont) {
      if (cont->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(cont);
        coefficient_construct_copy(ctx, cont, C);
      } else {
        integer_assign(ctx->K, &cont->value.num, &C->value.num);
      }
    }
    if (pp) {
      if (pp->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(pp);
        coefficient_construct_from_int(ctx, pp, 1);
      } else {
        integer_assign_int(ctx->K, &pp->value.num, 1);
      }
    }
    break;
  case COEFFICIENT_POLYNOMIAL:
  {
    int i;
    coefficient_t gcd;
    // Compute the gcd of coefficients starting with LC
    coefficient_construct_copy(ctx, &gcd, coefficient_lc(C));
    // Make if positive in case it's the only one
    if (coefficient_lc_sgn(ctx, &gcd) < 0) {
      coefficient_neg(ctx, &gcd, &gcd);
    }
    // Compute the rest of the gcd
    for (i = SIZE(C)-2; i >= 0 ; -- i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_gcd(ctx, &gcd, &gcd, COEFF(C, i));
        if (coefficient_is_one(ctx, &gcd)) {
          break;
        }
      }
    }
    // GCD is positive, so if the leading coefficient of C is negative, flip it
    if (coefficient_lc_sgn(ctx, C) < 0) {
      coefficient_neg(ctx, &gcd, &gcd);
    }

    if (pp) {
      // Now compute the pp
      coefficient_div(ctx, pp, C, &gcd);
      assert(coefficient_is_normalized(ctx, pp));
    }
    if (cont) {
      coefficient_swap(&gcd, cont);
      assert(coefficient_is_normalized(ctx, cont));
    }
    coefficient_destruct(&gcd);
    break;
  }
  default:
    assert(0);
    break;
  }

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_pp_cont() => ");
    if (pp) { tracef("pp = "); coefficient_print(ctx, pp, trace_out); tracef("\n"); }
    if (cont) { tracef("cont = "); coefficient_print(ctx, cont, trace_out); tracef("\n"); }
  }
}

void coefficient_cont(const lp_polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C) {
  coefficient_pp_cont(ctx, 0, cont, C);
}

void coefficient_pp(const lp_polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C) {
  coefficient_pp_cont(ctx, pp, 0, C);
}

lp_polynomial_vector_t* coefficient_mgcd(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2, const lp_assignment_t* m) {

  // Only for polynomials of the same type
  assert(C1->type == COEFFICIENT_POLYNOMIAL);
  assert(C2->type == COEFFICIENT_POLYNOMIAL);
  assert(coefficient_top_variable(C1) == coefficient_top_variable(C2));

  TRACE("coefficient", "coefficient_lcm()\n");

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
  }

  lp_variable_t x = coefficient_top_variable(C1);

  coefficient_t A, B, P, R, cont;
  coefficient_construct_copy(ctx, &A, C1);
  coefficient_construct_copy(ctx, &B, C2);
  coefficient_construct(ctx, &P);
  coefficient_construct(ctx, &R);
  coefficient_construct(ctx, &cont);

  lp_polynomial_vector_t* assumptions = lp_polynomial_vector_new(ctx);

  // Get the reductums of A and B
  coefficient_reductum_m(ctx, &A, &A, m, assumptions);
  coefficient_reductum_m(ctx, &B, &B, m, assumptions);

  // Get the primitive parts (reductum includes the sign of cont)
  coefficient_pp_cont(ctx, &A, &cont, &A);
  lp_polynomial_vector_push_back_coeff_prime(assumptions, &cont);
  coefficient_pp_cont(ctx, &B, &cont, &B);
  lp_polynomial_vector_push_back_coeff_prime(assumptions, &cont);

  // If one of the coefficient reduces to a constant, we're done
  if (coefficient_top_variable(&A) != x || coefficient_top_variable(&B) != x) {
    return assumptions;
  }

  // Swap A and B if def(A) < deg(B)
  if (coefficient_degree(&A) < coefficient_degree(&B)) {
    coefficient_swap(&A, &B);
  }

  //
  // We compute the reduction of A and B in Z[y, x], i.e.
  //
  //   P*A = Q*B + R
  //
  // with P in Z[y], Q in Z[y, x], and deg(R) < deg(B) or deg(R) == 0.
  //
  // We keep the accumulating the assumptions of the reduction and keep A, B, R
  // such reduced my model and primitive.
  //
  do {

    if (trace_is_enabled("coefficient::mgcd")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
    }


    // One step reduction, we get P*A = Q*B + R
    // If A, B have a common zero, this is also a zero of R (if R is in x)
    // If B, R have a common zero, this is also a zero of A if P != 0
    coefficient_reduce(ctx, &A, &B, &P, 0, &R, REMAINDERING_LCM_SPARSE);
//    if (!coefficient_is_constant(&P)) {
//      lp_polynomial_vector_push_back_coeff_prime(assumptions, &P);
//    }

    // Reduce R and pp
    if (!coefficient_is_constant(&R)) {
      coefficient_reductum_m(ctx, &R, &R, m, assumptions);
    }
    coefficient_pp_cont(ctx, &R, &cont, &R);
    lp_polynomial_vector_push_back_coeff_prime(assumptions, &cont);

    // We continue if we didn't get a constant
    int cmp_type = coefficient_cmp_type(ctx, &B, &R);
    if (cmp_type == 0) {
       // A = B, B = R (already reduced)
      coefficient_swap(&A, &B);
      coefficient_swap(&B, &R);
    } else {
      // Got to the GCD, but we need to maintain the sign of R
      if (!coefficient_is_constant(&R)) {
        lp_polynomial_vector_push_back_coeff_prime(assumptions, &R);
      }
      break;
    }
  } while (1);

  // Return the assumptions
  return assumptions;
}
