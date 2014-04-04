/*
 * output.c
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#include "utils/output.h"

#include <assert.h>

poly_output_language_t output_language = OUTPUT_LATEX;

void set_output_language(poly_output_language_t lang) {
  output_language = lang;
}

const char* get_power_symbol(void) {
 switch (output_language) {
 case OUTPUT_LATEX:
   return "^";
   break;
 case OUTPUT_SMT2:
   assert(0);
   return "^";
   break;
 case OUTPUT_PYTHON:
   return "**";
   break;
 default:
   assert(0);
   break;
 }
 return "^";
}



