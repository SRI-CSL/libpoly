/*
 * univariate_polynomial_factorization.h
 *
 *  Created on: Nov 8, 2013
 *      Author: dejan
 */

#pragma once

#include <stdio.h>

#include "upolynomial/upolynomial.h"
#include "upolynomial/factors.h"

/**
 * Factors the given polynomial into square-free factors. Polynomial f should be
 * monic if in Z_p, or primitive if in Z.
 */
lp_upolynomial_factors_t* upolynomial_factor_square_free(const lp_upolynomial_t* f);

/**
 * Factors the given polynomial into a distinct degree factorization.
 * Polynomial f should in Z_p, square-free, and monic.
 *
 * The output of the function is a factorization where each factor is associated
 * with a distinct degree of each of its sub-factors.
 */
lp_upolynomial_factors_t* upolynomial_factor_distinct_degree(const lp_upolynomial_t* f);

/**
 * Factors the given polynomial using the algorithm of Berlekamp. The algorithm
 * assumes that p is in a ring Z_p for some prime p.
 */
lp_upolynomial_factors_t* upolynomial_factor_Zp(const lp_upolynomial_t* f);

/**
 * Factors the given polynomial using Hansel lifting.
 */
lp_upolynomial_factors_t* upolynomial_factor_Z(const lp_upolynomial_t* f);
