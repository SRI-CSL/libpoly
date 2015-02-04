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

//
// Definitions of all relevant types
//

typedef struct lp_upolynomial_struct lp_upolynomial_t;
typedef struct lp_upolynomial_factors_struct lp_upolynomial_factors_t;
typedef struct lp_algebraic_number_struct lp_algebraic_number_t;
typedef struct lp_polynomial_context_struct lp_polynomial_context_t;
typedef struct lp_polynomial_struct lp_polynomial_t;
typedef struct lp_assignment_struct lp_assignment_t;

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



