/*
 * factorization.c
 *
 *  Created on: Mar 26, 2014
 *      Author: dejan
 */

#include "polynomial/factorization.h"

#include "polynomial/gcd.h"
#include "polynomial/output.h"

#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <malloc.h>

STAT_DECLARE(int, coefficient, factor_square_free)
STAT_DECLARE(int, coefficient, factor_square_free_pp)

#define INITIAL_FACTORS_CAPACITY 5

void coefficient_factors_construct(coefficient_factors_t* factors) {
  factors->size = 0;
  factors->capacity = INITIAL_FACTORS_CAPACITY;
  factors->factors = malloc(sizeof(coefficient_t)*INITIAL_FACTORS_CAPACITY);
  factors->multiplicities = malloc(sizeof(size_t)*INITIAL_FACTORS_CAPACITY);
}

void coefficient_factors_destruct(coefficient_factors_t* factors) {
  if (factors->factors) {
    size_t i;
    for (i = 0; i < factors->size; ++ i) {
      coefficient_destruct(factors->factors + i);
    }
    free(factors->factors);
    free(factors->multiplicities);
  }
}

void coefficient_factors_add(const lp_polynomial_context_t* ctx, coefficient_factors_t* factors, const coefficient_t* C, size_t multiplicity) {
  if (factors->size == factors->capacity) {
    factors->capacity *= 2;
    factors->factors = realloc(factors->factors, sizeof(coefficient_t)*factors->capacity);
    factors->multiplicities = realloc(factors->multiplicities, sizeof(size_t)*factors->capacity);
  }

  factors->multiplicities[factors->size] = multiplicity;
  coefficient_construct_copy(ctx, factors->factors + factors->size, C);
  factors->size ++;
}

int coefficient_factors_print(const lp_polynomial_context_t* ctx, const coefficient_factors_t* factors, FILE* out) {
  int ret = 0;

  fprintf(out, "[");

  size_t i;
  for (i = 0; i < factors->size; ++ i) {
    if (i) {
      ret += fprintf(out, ", ");
    }
    ret += fprintf(out, "(");
    ret += coefficient_print(ctx, factors->factors + i, out);
    ret += fprintf(out, ", %zu)", factors->multiplicities[i]);
  }

  fprintf(out, "]");

  return ret;
}

/**
 * We are given a polynomial f and we will return it's square-free factorization
 *
 *  f = \prod_{k} f_k^k
 *
 * where each f_k is square free and all gcd(f_i, f_j) = 1.
 *
 * First, note that if f has a square factor, then the gcd(f, f') != 1. This is
 * because if f = u^2*v, then f' = 2*u*u'*v + u^2*v, so gcd(f, f') would include
 * at least u != 1.
 *
 * Now,consider f'
 *
 *  f' = \sum_{k} k*f_k^{k-1}*f_k'*\prod_{i!=k} f_i^i (if we are in Z_p)
 *     = \sum_{p \ndiv k} k*f_k^{k-1}*f_k'*\prod_{i!=k}
 *
 * It is easy to see that f' is divisible by each f_k^{k-1} for k such that
 * p \ndiv k, but not f_k^k. It is also divisible by f_k^k for k such that p
 * \div k. Since these are the only possible factors of f, then
 *
 *  P = gcd(f, f') = \prod_{p \ndiv k} f_k^{k-1} * \prod_{p \div k} f_k^k
 *
 * With P as above, then
 *
 *  L = f/P = \prod_{p \ndiv k} f_k               -> Linear product (can be 1)
 *
 * To compute the first factor we then compute (loop start)
 *
 *  R = gcd(P, L) = \prod_{k > 1 & p \ndiv k} f_k -> Rest of linear product
 *  O = L / R                                     -> First factor (output)
 *
 * To continue the loop, we set
 *
 *  P = P / R     -> Reduced powers of f_k with p \ndiv k, k > 1
 *  L = R         -> Linear product of f_k with p \ndiv k, k > 1
 *
 * And go on with the loop.
 *
 * The loop ends when L = 1 and then we know that P = \prod_{p \div k} f_k^k
 * which we know how to special case (if P != 1).
 */
void coefficient_factor_square_free_pp(const lp_polynomial_context_t* ctx, const coefficient_t* C, coefficient_factors_t* factors) {

  STAT(coefficient, factor_square_free_pp) ++;

  if (trace_is_enabled("factorization")) {
    tracef("coefficient_factor_square_free_pp("); coefficient_print(ctx, C, trace_out); tracef(")\n");
  }

  assert(C->type == COEFFICIENT_POLYNOMIAL);

  // Derivative
  coefficient_t C_d;
  coefficient_construct(ctx, &C_d);
  coefficient_derivative(ctx, &C_d, C);

  if (coefficient_is_zero(ctx, &C_d)) {
      assert(ctx->K && ctx->K->is_prime);
      // f' is zero for a non-zero polynomial => f has to be of the form
      // f = \sum a_k x^(p*d_k) = f_p(x^p) where f_p = \sum a_k x^d_k
      // we factor f_p and then return f_p(x^p)=(f_p)^p
      size_t p = (size_t) integer_to_int(&ctx->K->M);
      coefficient_t C_p;
      coefficient_construct_copy(ctx, &C_p, C);
      coefficient_div_degrees(ctx, &C_p, p);
      size_t i = factors->size;
      coefficient_factor_square_free_pp(ctx, &C_p, factors);
      // Adjust the multiplicities
      for (;i < factors->size; ++ i) {
        factors->multiplicities[i] *= p;
      }
      coefficient_destruct(&C_p);
    } else {

      // Degree of the factor
      size_t k = 1;

      // P = GCD(f, f')
      coefficient_t P;
      coefficient_construct(ctx, &P);
      coefficient_gcd(ctx, &P, C, &C_d);
      if (trace_is_enabled("factorization")) {
        tracef("P = "); coefficient_print(ctx, &P, trace_out); tracef("\n");
      }
      // L = f/P
      coefficient_t L;
      coefficient_construct(ctx, &L);
      coefficient_div(ctx, &L, C, &P);
      if (trace_is_enabled("factorization")) {
        tracef("L = "); coefficient_print(ctx, &L, trace_out); tracef("\n");
      }

      coefficient_t O, R;
      coefficient_construct(ctx, &R);
      coefficient_construct(ctx, &O);
      while (coefficient_degree(&L) > 0) {
        // R = gcd(P, L)
        coefficient_gcd(ctx, &R, &P, &L);
        if (trace_is_enabled("factorization")) {
          tracef("R = "); coefficient_print(ctx, &R, trace_out); tracef("\n");
        }
        // O = L / R (it can be constant if there is no factor of power k)
        if (coefficient_cmp(ctx, &L, &R)) {
          coefficient_div(ctx, &O, &L, &R);
          if (trace_is_enabled("factorization")) {
            tracef("O = "); coefficient_print(ctx, &O, trace_out); tracef("\n");
          }
          // Record the output
          coefficient_factors_add(ctx, factors, &O, k);
        }
        // P = P / R
        coefficient_div(ctx, &P, &P, &R);
        if (trace_is_enabled("factorization")) {
          tracef("P = "); coefficient_print(ctx, &P, trace_out); tracef("\n");
        }
        // L = R
        coefficient_swap(&L, &R);
        if (trace_is_enabled("factorization")) {
          tracef("L = "); coefficient_print(ctx, &L, trace_out); tracef("\n");
        }
        // Next degree
        k = k + 1;
      }

      // If P has content, it is a power of p
      if (coefficient_degree(&P) > 0) {
        size_t p = integer_to_int(&ctx->K->M);
        coefficient_div_degrees(ctx, &P, p);
        size_t i = factors->size;
        coefficient_factor_square_free_pp(ctx, &P, factors);
        for (; i < factors->size; ++ i) {
          factors->multiplicities[i] *= p;
        }
      }

      coefficient_destruct(&P);
      coefficient_destruct(&L);
      coefficient_destruct(&O);
      coefficient_destruct(&R);
    }


  coefficient_destruct(&C_d);

  if (trace_is_enabled("factorization")) {
    tracef("upolynomial_factor_square_free("); coefficient_print(ctx, C, trace_out); tracef(") = ");
    coefficient_factors_print(ctx, factors, trace_out); tracef("\n");
  }
}

void coefficient_factor_square_free(const lp_polynomial_context_t* ctx, const coefficient_t* C, coefficient_factors_t* factors) {

  STAT(coefficient, factor_square_free) ++;

  if (trace_is_enabled("factorization")) {
    tracef("coefficient_factor_square_free("); coefficient_print(ctx, C, trace_out); tracef(")\n");
  }

  coefficient_t C_pp, C_cont;
  coefficient_construct(ctx, &C_pp);
  coefficient_construct(ctx, &C_cont);

  // Get the content and primitive part
  coefficient_pp_cont(ctx, &C_pp, &C_cont, C);

  // Factor the content if not trivial
  if (coefficient_is_constant(&C_cont)) {
    if (!coefficient_is_one(ctx, &C_cont)) {
      coefficient_factors_add(ctx, factors, &C_cont, 1);
    }
  } else {
    coefficient_factor_square_free(ctx, &C_cont, factors);
  }

  // Factor the primitive part
  if (!coefficient_is_constant(&C_pp)) {
    coefficient_factor_square_free_pp(ctx, &C_pp, factors);
  }

  coefficient_destruct(&C_pp);
  coefficient_destruct(&C_cont);

  if (trace_is_enabled("factorization")) {
    tracef("coefficient_factor_square_free("); coefficient_print(ctx, C, trace_out); tracef(") =>");
    coefficient_factors_print(ctx, factors, trace_out); tracef("\n");

  }
}


