/*
 * debug_trace.h
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#pragma once

#include <stdio.h>

/** Interface to the tracing functionality */
typedef struct {

  void (*enable) (const char* tag);
  void (*disable) (const char* tag);
  int (*is_enabled) (const char* tag);

  void (*set_output) (FILE* file);

} debug_trace_ops_t;

/** Implementation of  the tracing functionality */
extern const debug_trace_ops_t debug_trace_ops;

