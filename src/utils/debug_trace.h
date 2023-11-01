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

#include <stdio.h>

/** Where the output goes (defaults to stderr). Use the macro below */
extern
FILE* trace_out_real;

#define trace_out (trace_out_real ? trace_out_real : stderr)

/** Print to the debug trace printf style */
#define tracef(...) fprintf(trace_out, __VA_ARGS__)

/** Set the output file for tracing */
void trace_set_output(FILE* file);

void trace_enable(const char* tag);

void trace_disable(const char* tag);

#ifndef NDEBUG

int trace_is_enabled(const char* tag);

#define TRACE(tag, ...) do { \
  if (trace_is_enabled(tag)) { \
    tracef(__VA_ARGS__); \
  } \
} while(0)

#define TRACE_CMD(tag, cmd) do { \
  if (trace_is_enabled(tag)) { \
    cmd; \
  } \
} while(0)

#else

static inline
int trace_is_enabled(const char* tag) {
  (void) tag;
  return 0;
}

#define TRACE(tag, ...) ((void)0)
#define TRACE_CMD(tag, cmd) ((void)0)
#endif
