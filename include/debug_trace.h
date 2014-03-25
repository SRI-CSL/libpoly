/*
 * debug_trace.h
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>

/**
 * Interface to the tracing functionality. The functionality is eabled only
 * in debug builds.
 */
typedef struct {

  /** Enable a given tag for tracing */
  void (*enable) (const char* tag);

  /** Disable the given that from tracing */
  void (*disable) (const char* tag);

  /** Check if the tag is enabled for tracing */
  int (*is_enabled) (const char* tag);

  /** Set the output trace file (defaults to stderr) */
  void (*set_output) (FILE* file);

} debug_trace_ops_t;

/** Implementation of  the tracing functionality */
extern const debug_trace_ops_t debug_trace_ops;

