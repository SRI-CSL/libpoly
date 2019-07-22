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

#include "utils/output.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

__thread
lp_output_language_t output_language = LP_OUTPUT_LATEX;

__thread
char* upolynomial_var_symbol = 0;

void set_output_language(lp_output_language_t lang) {
  output_language = lang;
}

void set_upolynomial_var_symbol(const char* x) {
  if (upolynomial_var_symbol) {
    free(upolynomial_var_symbol);
  }
  upolynomial_var_symbol = strdup(x);
}

const char* get_upolynomial_var_symbol() {
  if (upolynomial_var_symbol) {
    return upolynomial_var_symbol;
  } else {
    return "x";
  }
}

const char* get_power_symbol(void) {
 switch (output_language) {
 case LP_OUTPUT_LATEX:
   return "^";
   break;
 case LP_OUTPUT_MATHEMATICA:
   return "^";
   break;
 case LP_OUTPUT_SMT2:
   assert(0);
   return "^";
   break;
 case LP_OUTPUT_PYTHON:
   return "**";
   break;
 default:
   assert(0);
   break;
 }
 return "^";
}
