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

#include "utils/debug_trace.h"

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

__thread
char* tags_to_trace[1000];

__thread
size_t tags_to_trace_size = 0;

void trace_enable(const char* tag) {
#ifndef NDEBUG
  tags_to_trace[tags_to_trace_size++] = strdup(tag);
#else
  (void) tag;
#endif
}

#ifndef NDEBUG
int trace_is_enabled(const char* tag) {
  unsigned i;
  for (i = 0; i < tags_to_trace_size; ++ i) {
    if (strcmp(tag, tags_to_trace[i]) == 0) {
      return i+1;
    }
  }
  return 0;
}
#endif

void trace_disable(const char* tag) {
#ifndef NDEBUG
  int i = trace_is_enabled(tag) - 1;
  if (i >= 0) {
    free(tags_to_trace[i]);
    tags_to_trace[i] = tags_to_trace[--tags_to_trace_size];
  }
#else
  (void) tag;
#endif
}

__thread
FILE* trace_out_real = 0;

void trace_set_output(FILE* file) {
  trace_out_real = file;
}
