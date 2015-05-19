/*
 * poly.h
 *
 *  Created on: Mar 24, 2014
 *      Author: dejan
 */

#pragma once

#define LIBPOLY_VERSION_MAJOR 0
#define LIBPOLY_VERSION_MINOR 0
#define LIBPOLY_VERSION_GIT 0

#include "output_language.h"

#include <stdio.h>
#include <gmp.h>

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

/**
 * Function to initialize the library (statistics, logging). Should be called
 * before any other function.
 */
void lp_init();

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



