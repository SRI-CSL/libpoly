/*
 * upolynomial_factors.h
 *
 *  Created on: Feb 4, 2015
 *      Author: dejan
 */

#pragma once

#include "poly.h"
#include "integer.h"

/** Construct the factors with no factors, and constant 1. */
lp_upolynomial_factors_t* lp_upolynomial_factors_construct(void);

/**
 * Free the memory and the factor polynomials. If destruct_factors is true
 * then the individual factors are also destructed (should be, unless you
 * copied the factors somewhere else).
 */
void lp_upolynomial_factors_destruct(lp_upolynomial_factors_t* f, int destruct_factors);

/** Clear the factor polynomials */
void lp_upolynomial_factors_clear(lp_upolynomial_factors_t* f);

/** Swap the two factorizations */
void lp_upolynomial_factors_swap(lp_upolynomial_factors_t* f1, lp_upolynomial_factors_t* f2);

/** Get the number of factors */
size_t lp_upolynomial_factors_size(const lp_upolynomial_factors_t* f);

/** Get a factor with the given index i < size() */
lp_upolynomial_t* lp_upolynomial_factors_get_factor(lp_upolynomial_factors_t* f, size_t i, size_t* multiplicity);

/** Returns the constant of the factorization */
const lp_integer_t* lp_upolynomial_factors_get_constant(const lp_upolynomial_factors_t* f);

/** Add a factor with the given degree */
void lp_upolynomial_factors_add(lp_upolynomial_factors_t* f, lp_upolynomial_t* p, size_t d);

/** Print the factors */
int lp_upolynomial_factors_print(const lp_upolynomial_factors_t* f, FILE* out);

/** Get the ring */
lp_int_ring lp_upolynomial_factors_ring(const lp_upolynomial_factors_t* f);

/**
 * Set the ring of all polynomials to K. This is only possible if K is
 * "larger" than the existing ring.
 */
void lp_upolynomial_factors_set_ring(lp_upolynomial_factors_t* f, lp_int_ring K);


