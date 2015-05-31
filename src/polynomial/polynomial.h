/*
 * polynomial_internal.h
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#pragma once

#include <polynomial.h>

#include "polynomial/coefficient.h"

struct lp_polynomial_struct {
  /** The actual polynomial representation (so we can use it as a coefficient) */
  coefficient_t data;
  /** Hash */
  size_t hash;
  /** Is this an external polynomial (needs checks on function entry) */
  char external;
  /** Context of the polynomial */
  const lp_polynomial_context_t* ctx;
};

/** Construct from coefficient */
void lp_polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from);

