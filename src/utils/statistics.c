/*
 * statistics.c
 *
 *  Created on: Dec 4, 2013
 *      Author: dejan
 */

#include "statistics.h"

#include <string.h>
#include <stdlib.h>

#define INITIAL_STATS 100

typedef struct int_stats_struct {
  size_t count;
  int values[INITIAL_STATS];
  char* names[INITIAL_STATS];
} int_stats_t;

/** Integer statistics */
static int_stats_t int_stats;

int* stats_register_int(const char* name) {
  size_t i = int_stats.count ++;
  int_stats.values[i] = 0;
  int_stats.names[i] = strdup(name);
  return int_stats.values + i;
}

void stats_print(FILE* out) {
  size_t i;
  for (i = 0; i < int_stats.count; ++ i) {
    fprintf(out, "%s = %d\n", int_stats.names[i], int_stats.values[i]);
  }
}

void stats_construct(void) {
  int_stats.count = 0;
}

void stats_destruct(void) {
  size_t i;
  for (i = 0;i < int_stats.count; ++ i) {
    free(int_stats.names[i]);
  }
}

