/*
 * upolynomial.h
 *
 *  Created on: Feb 4, 2015
 *      Author: dejan
 */

#pragma once

#include <integer.h>
#include "upolynomial/umonomial.h"

/**
 * A polynomial is the ring, number of monomials and the monomials.
 */
struct lp_upolynomial_struct {
  /** The ring of coefficients */
  lp_int_ring K;
  /** The number of monomials */
  size_t size;
  /** The monomials */
  umonomial_t monomials[];
};
