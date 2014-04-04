/*
 * output_language.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

enum poly_output_language_enum {
  /** Output in latex format, i.e. x^2 - 1 */
  OUTPUT_LATEX,
  /** Output in SMT2 format, i.e. (+ (* x x) (- 1)) */
  OUTPUT_SMT2,
  /** Output in Python format, i.e. x**2 - 1 */
  OUTPUT_PYTHON,
};

typedef enum poly_output_language_enum poly_output_language_t;
