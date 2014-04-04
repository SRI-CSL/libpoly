/*
 * debug_trace_internal.h
 *
 *  Created on: Mar 21, 2014
 *      Author: dejan
 */

#pragma once

#include <stdio.h>

/** Where the output goes (defaults to */
extern FILE* trace_out;

/** Print to the debug trace printf style */
#define tracef(...) fprintf(trace_out, __VA_ARGS__);

#ifndef NDEBUG

void trace_enable(const char* tag);

void trace_disable(const char* tag);

int trace_is_enabled(const char* tag);

void trace_set_output(FILE* file);

#define TRACE(tag, ...) { \
  if (trace_is_enabled(tag)) { \
    tracef(__VA_ARGS__); \
  } \
}

#define TRACE_CMD(tag, cmd) { \
  if (trace_is_enabled(tag)) { \
    cmd; \
  } \
} \

#else

static inline
int trace_is_enabled(const char* tag) {
  (void) tag;
  return 0;
}

#define TRACE(tag, fmt, ...)
#define TRACE_CMD(tag, cmd)
#endif
