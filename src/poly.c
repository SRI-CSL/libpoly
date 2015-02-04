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

const lp_ops_t lp_poly_ops = {
    trace_enable,
    trace_disable,
    trace_set_output,
    stats_print,
    set_output_language
};
