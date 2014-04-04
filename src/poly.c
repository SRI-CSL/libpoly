/*
 * poly.c
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#include <poly.h>

#include "utils/statistics.h"
#include "utils/debug_trace.h"

void poly_init(void) {
  stats_construct();
  trace_set_output(stderr);
}

void poly_set_output_language(poly_output_language_t lang) {
  (void) lang;
}

const poly_ops_t poly_ops = {
    poly_init,
    trace_enable,
    trace_disable,
    trace_set_output,
    stats_print,
    poly_set_output_language
};
