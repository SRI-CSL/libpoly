/*
 * factors.h
 *
 *  Created on: Feb 4, 2015
 *      Author: dejan
 */

#pragma once

#include <integer.h>
#include <upolynomial.h>

struct lp_upolynomial_factors_struct {
  /** Constant factor */
  lp_integer_t constant;
  /** Number of actual factors */
  size_t size;
  /** Size of the factors array */
  size_t capacity;
  /** The irreducible factors */
  lp_upolynomial_t** factors;
  /** The multiplicity of individual factors */
  size_t* multiplicities;
};
