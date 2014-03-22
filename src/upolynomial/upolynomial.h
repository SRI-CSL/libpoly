/*
 * univariate_polynomial_internal.h
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#pragma once

#include <upolynomial.h>

#include "umonomial.h"

/**
 * A polynomial is the ring, number of monomials and the monomials.
 */
typedef struct upolynomial_struct {

  /** The ring of coefficients */
  int_ring K;
  /** The number of monomials */
  size_t size;
  /** The monomials */
  umonomial_t monomials[0];

} upolynomial_t;


typedef struct upolynomial_factors_struct {

  /** Constant factor */
  integer_t constant;
  /** Number of actual factors */
  size_t size;
  /** Size of the factors array */
  size_t capacity;
  /** The irreducible factors */
  upolynomial_t** factors;
  /** The multiplicity of individual factors */
  size_t* multiplicities;

} upolynomial_factors_t;


