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
// Definitions of all relevan types
//

typedef struct upolynomial_struct upolynomial_t;
typedef struct upolynomial_factors_struct upolynomial_factors_t;
typedef struct algebraic_number_struct algebraic_number_t;
typedef struct polynomial_context_struct polynomial_context_t;
typedef struct polynomial_struct polynomial_t;
typedef struct assignment_struct assignment_t;

typedef struct {

  /** Initialize the library (should be called before anything else) */
  void (*init) (void);

  /** Enable a given tag for tracing */
  void (*trace_enable) (const char* tag);

  /** Disable the given that from tracing */
  void (*trace_disable) (const char* tag);

  /** Set the output trace file (defaults to stderr) */
  void (*trace_set_output) (FILE* file);

  /** Print statistics to a file */
  void (*stats_print) (FILE* file);

  /** Set the output language */
  void (*set_output_language) (poly_output_language_t lang);

} poly_ops_t;

extern const poly_ops_t poly_ops;
