/*
 * debug_trace_internal.h
 *
 *  Created on: Mar 21, 2014
 *      Author: dejan
 */

#pragma once

#include <debug_trace.h>

/** Where the output goes */
extern FILE* trace_out;

/** Print to the debug trace printf style */
#define tracef(...) { \
  if (trace_out) { \
    fprintf(trace_out, __VA_ARGS__); \
  } else { \
    fprintf(stderr, __VA_ARGS__); \
  } \
}

#ifndef NDEBUG

#define TRACE(tag, ...) { \
  if (debug_trace_ops.is_enabled(tag)) { \
    tracef(__VA_ARGS__); \
  } \
}

#define TRACE_CMD(tag, cmd) { \
  if (debug_trace_ops.is_enabled(tag)) { \
    cmd; \
  } \
} \

#else
#define TRACE(tag, fmt, ...)
#define TRACE_CMD(tag, cmd)
#endif
