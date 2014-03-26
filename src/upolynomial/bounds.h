/*
 * univariate_polynomial_bounds.h
 *
 *  Created on: Dec 6, 2013
 *      Author: dejan
 */

#pragma once

#include "number/integer.h"
#include "upolynomial/upolynomial.h"

/**
 * Upper bound on the modulus of the roots of f.
 *
 * Given
 *
 *  f = a_n x^n + ... + a_0, with a_n != 0
 *
 * the bound is
 *
 *  B = 1 + max(a_0, ..., a_{n-1})/|a_n|
 */
void upolynomial_root_bound_cauchy(const upolynomial_t* f, integer_t* B);

/**
 * Computes the bound on the size of coefficients of any polynomial g with
 * deg(g) <= n that divides f.
 *
 * Given
 *
 *  f = a_n x^n + ... + a_0, with a_n != 0
 *
 * the bound is
 *
 *  B = 2^d*norm(f)
 */
void upolynomial_factor_bound_landau_mignotte(const upolynomial_t* f, size_t n, integer_t* B);
