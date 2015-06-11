/*
 * poly.c
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#include <poly.h>

#include "utils/statistics.h"
#include "utils/debug_trace.h"
#include "utils/output.h"


void lp_init(void) {
  stats_construct();
  trace_set_output(stderr);
}

void lp_trace_enable(const char* tag) {
  trace_enable(tag);
}

void lp_trace_disable(const char* tag) {
  trace_disable(tag);
}

void lp_trace_set_output(FILE* file) {
  trace_set_output(file);
}

void lp_stats_print(FILE* file) {
  stats_print(file);
}

void lp_set_output_language(lp_output_language_t lang) {
  set_output_language(lang);
}

void lp_set_upolynomial_var_symbol(const char* x) {
  set_upolynomial_var_symbol(x);
}
