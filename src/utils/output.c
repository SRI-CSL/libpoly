/*
 * output.c
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#include "utils/output.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

lp_output_language_t output_language = LP_OUTPUT_LATEX;

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



