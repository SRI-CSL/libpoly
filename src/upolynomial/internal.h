/*
 * univariate_polynomial_internal.h
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#pragma once

#include "number/integer.h"

/**
 * A monomial is a coefficient and the degree.
 */
typedef struct umonomial_struct {

  /** Degree of the monomial */
  size_t degree;
  /** Coefficient with the monomial */
  integer_t coefficient;

} umonomial_t;

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

void umonomial_construct(int_ring K, umonomial_t* m, size_t degree, const integer_t* coefficient);

void umonomial_construct_from_int(int_ring K, umonomial_t* m, size_t degree, long coefficient);

void umonomial_construct_copy(int_ring K, umonomial_t* m, const umonomial_t* from);

void umonomial_destruct(umonomial_t* m);

int umonomial_print(const umonomial_t* m, FILE* out);
