/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "version.h"
#include "output_language.h"

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

//
// Definitions of all relevant types
//

/** Rationals are GMP rationals */
typedef __mpq_struct lp_rational_t;

/** Integers are GMP integers */
typedef __mpz_struct lp_integer_t;


typedef size_t lp_variable_t;

#define lp_variable_null ((lp_variable_t)(-1))

typedef struct lp_variable_db_struct lp_variable_db_t;
typedef struct lp_variable_list_struct lp_variable_list_t;
typedef struct lp_variable_order_struct lp_variable_order_t;
typedef struct lp_variable_order_ops_struct lp_variable_order_ops_t;

typedef struct lp_upolynomial_struct lp_upolynomial_t;
typedef struct lp_upolynomial_factors_struct lp_upolynomial_factors_t;

typedef struct lp_polynomial_context_struct lp_polynomial_context_t;
typedef struct lp_polynomial_struct lp_polynomial_t;

typedef struct lp_algebraic_number_struct lp_algebraic_number_t;
typedef struct lp_value_struct lp_value_t;
typedef struct lp_assignment_struct lp_assignment_t;

typedef struct lp_rational_interval_struct lp_rational_interval_t;
typedef struct lp_dyadic_interval_struct lp_dyadic_interval_t;
typedef struct lp_interval_struct lp_interval_t;

typedef struct lp_feasibility_set_struct lp_feasibility_set_t;
typedef struct lp_polynomial_hash_set_struct lp_polynomial_hash_set_t;
typedef struct lp_polynomial_vector_struct lp_polynomial_vector_t;

/** Enable a given tag for tracing */
void lp_trace_enable(const char* tag);

/** Disable the given that from tracing */
void lp_trace_disable(const char* tag);

/** Set the output trace file (defaults to stderr) */
void lp_trace_set_output(FILE* file);

/** Print statistics to a file */
void lp_stats_print(FILE* file);

/** Set the output language */
void lp_set_output_language(lp_output_language_t lang);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
