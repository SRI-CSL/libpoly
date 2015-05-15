/*
 * polynomial_context.h
 *
 *  Created on: May 15, 2015
 *      Author: dejan
 */

#pragma once

#include <poly.h>

/** Get a temp variable */
lp_variable_t lp_polynomial_context_get_temp_variable(const lp_polynomial_context_t* ctx_const);

/** Release the variable (has to be the last one obtained and not released */
void lp_polynomial_context_release_temp_variable(const lp_polynomial_context_t* ctx_const, lp_variable_t x);
