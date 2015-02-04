/*
 * factorization.h
 *
 *  Created on: Mar 26, 2014
 *      Author: dejan
 */

#pragma once

#include "polynomial/coefficient.h"

#include <stdio.h>

struct coefficient_factors_struct {
  size_t size;
  size_t capacity;
  coefficient_t* factors;
  size_t* multiplicities;
};

typedef struct coefficient_factors_struct coefficient_factors_t;

void coefficient_factors_construct(coefficient_factors_t* factors);
void coefficient_factors_destruct(coefficient_factors_t* factors);
void coefficient_factors_add(const lp_polynomial_context_t* ctx, coefficient_factors_t* factors, const coefficient_t* C, size_t multiplicity);
int coefficient_factors_print(const lp_polynomial_context_t* ctx, const coefficient_factors_t* factors, FILE* out);

/** Factors the given coefficient into square-free factors. */
void coefficient_factor_square_free(const lp_polynomial_context_t* ctx, const coefficient_t* C, coefficient_factors_t* factors);

