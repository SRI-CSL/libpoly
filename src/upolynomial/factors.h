/*
 * factors.h
 *
 *  Created on: Mar 25, 2014
 *      Author: dejan
 */

#pragma once

#include <upolynomial.h>

#include "number/integer.h"

struct upolynomial_factors_struct {
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
};

upolynomial_factors_t* upolynomial_factors_construct(void);

void upolynomial_factors_swap(upolynomial_factors_t* f1, upolynomial_factors_t* f2);

void upolynomial_factors_clear(upolynomial_factors_t* f);

void upolynomial_factors_destruct(upolynomial_factors_t* f, int destruct_factors);

size_t upolynomial_factors_size(const upolynomial_factors_t* f);

upolynomial_t* upolynomial_factors_get_factor(upolynomial_factors_t* f, size_t i, size_t* d);

const integer_t* upolynomial_factors_get_constant(const upolynomial_factors_t* f);

void upolynomial_factors_add(upolynomial_factors_t* f, upolynomial_t* p, size_t d);

int upolynomial_factors_print(const upolynomial_factors_t* f, FILE* out);

int_ring upolynomial_factors_ring(const upolynomial_factors_t* f);

void upolynomial_factors_set_ring(upolynomial_factors_t* f, int_ring K);
