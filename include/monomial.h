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

#pragma once

#include "integer.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Power of single variable.
 */
typedef struct {
  lp_variable_t x;
  size_t d;
} power_t;

/**
 * A monomial of the form a*x_1^d_1*...*x_size*d_size.
 */
typedef struct {
  /** The coefficient */
  lp_integer_t a;
  /** Number of variables */
  size_t n;
  /** The capacity of the power array */
  size_t capacity;
  /** Array of variable powers */
  power_t* p;
} lp_monomial_t;

/** Construct an empty monomial (0) */
void lp_monomial_construct(const lp_polynomial_context_t* ctx, lp_monomial_t* m);

/** Construct an ordered copy of the monomial */
void lp_monomial_construct_copy(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_monomial_t* from, int sort);

/** Set the constant of the monomial */
void lp_monomial_set_coefficient(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_integer_t* a);

/** Destruct the monomial */
void lp_monomial_destruct(lp_monomial_t* m);

/** Clear the monomial to 0 */
void lp_monomial_clear(const lp_polynomial_context_t* ctx, lp_monomial_t* m);

/** Assign another monomial */
void lp_monomial_assign(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_monomial_t* from, int sort);

/** Add a variable power to the end of the monomial */
void lp_monomial_push(lp_monomial_t* m, lp_variable_t x, size_t d);

/** Remove a variable power from the end of the monomial */
void lp_monomial_pop(lp_monomial_t* m);

/** Get the gcd of two monomials */
void lp_monomial_gcd(const lp_polynomial_context_t* ctx, lp_monomial_t* gcd, const lp_monomial_t* m1, const lp_monomial_t* m2);

/** Print the monomial */
int lp_monomial_print(const lp_polynomial_context_t* ctx, const lp_monomial_t* m, FILE* out);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
