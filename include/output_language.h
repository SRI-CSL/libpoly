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

#ifdef __cplusplus
extern "C" {
#endif

enum lp_output_language_enum {
  /** Output in latex format, i.e. x^2 - 1 */
  LP_OUTPUT_LATEX,
  /** Output in SMT2 format, i.e. (+ (* x x) (- 1)) */
  LP_OUTPUT_SMT2,
  /** Output in Python format, i.e. x**2 - 1 */
  LP_OUTPUT_PYTHON,
  /** Output in Mathematica format */
  LP_OUTPUT_MATHEMATICA
};

typedef enum lp_output_language_enum lp_output_language_t;

#ifdef __cplusplus
} /* close extern "C" { */
#endif
