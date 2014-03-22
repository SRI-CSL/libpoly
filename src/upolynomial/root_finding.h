/*
 * univariate_polynomial_root_finding.h
 *
 *  Created on: Jan 13, 2014
 *      Author: dejan
 */

#pragma once

#include "upolynomial/upolynomial.h"
#include "upolynomial/internal.h"
#include "upolynomial/upolynomial_dense.h"

#include "number/dyadic_rational.h"
#include "number/algebraic_number.h"

#include "interval/interval.h"

/**
 * Compute the Sturm sequence of f. A Sturm sequence is a sequence
 *
 *  S[0] = f
 *  S[1] = f'
 *
 *  a*S[i-2] = Q*S[i-1] + b*S[i]
 *
 * with a*b < 0.
 */
void upolynomial_compute_sturm_sequence(const upolynomial_t* f, upolynomial_dense_t* S, size_t* size);

/**
 * Count the number of real roots that the polynomial f has in the given open
 * interval. The polynomial f should be square-free.
 */
int upolynomial_roots_count_sturm(const upolynomial_t* f, const interval_t* interval);

/**
 * Isolate the root intervals of the polynomial f and construct the resulting
 * numbers into the given array (should be at least deg(f) size). The size
 * will be updated to the number of roots.
 */
void upolynomial_roots_isolate_sturm(const upolynomial_t* f, algebraic_number_t* roots, size_t* roots_size);
