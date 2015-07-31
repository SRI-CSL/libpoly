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

#include <poly.h>

#include "utils/statistics.h"
#include "utils/debug_trace.h"
#include "utils/output.h"

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
