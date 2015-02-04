/*
 * factors.h
 *
 *  Created on: Mar 25, 2014
 *      Author: dejan
 */

#pragma once

#include <upolynomial.h>

#include "number/integer.h"

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

lp_upolynomial_factors_t* upolynomial_factors_construct(void);

void upolynomial_factors_swap(lp_upolynomial_factors_t* f1, lp_upolynomial_factors_t* f2);

void upolynomial_factors_clear(lp_upolynomial_factors_t* f);

void upolynomial_factors_destruct(lp_upolynomial_factors_t* f, int destruct_factors);

size_t upolynomial_factors_size(const lp_upolynomial_factors_t* f);

lp_upolynomial_t* upolynomial_factors_get_factor(lp_upolynomial_factors_t* f, size_t i, size_t* d);

const lp_integer_t* upolynomial_factors_get_constant(const lp_upolynomial_factors_t* f);

void upolynomial_factors_add(lp_upolynomial_factors_t* f, lp_upolynomial_t* p, size_t d);

lp_int_ring upolynomial_factors_ring(const lp_upolynomial_factors_t* f);

void upolynomial_factors_set_ring(lp_upolynomial_factors_t* f, lp_int_ring K);
