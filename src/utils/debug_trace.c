/*
 * debug_trace.c
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#include "utils/debug_trace.h"

#include <string.h>
#include <malloc.h>
#include <stdarg.h>

char* tags_to_trace[1000];
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

FILE* trace_out = 0;

void trace_set_output(FILE* file) {
  trace_out = file;
}
