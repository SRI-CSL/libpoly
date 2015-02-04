/*
 * output.h
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#pragma once

#include "polynomial/coefficient.h"
#include "utils/output.h"

/** Print the monomial */
int monomial_print(const lp_polynomial_context_t* ctx, const monomial_t* m, FILE* out);

/** Prints the coefficient to the given stream. */
int coefficient_print(const lp_polynomial_context_t* ctx, const coefficient_t* C, FILE* out);

/** Returns the string representation of the coefficient. */
char* coefficient_to_string(const lp_polynomial_context_t* ctx, const coefficient_t* C);

