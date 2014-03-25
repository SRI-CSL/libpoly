/*
 * statistics.c
 *
 *  Created on: Dec 4, 2013
 *      Author: dejan
 */

#include "statistics.h"

#include <string.h>
#include <malloc.h>

#define INITIAL_STATS 100

typedef struct int_stats_struct {
  size_t count;
  int values[INITIAL_STATS];
  char* names[INITIAL_STATS];
} int_stats_t;

/** Integer statistics */
static int_stats_t int_stats;

static int* stats_register_int(const char* name) {
  size_t i = int_stats.count ++;
  int_stats.values[i] = 0;
  int_stats.names[i] = strdup(name);
  return int_stats.values + i;
}

__attribute__ (( __constructor__ (101) ))
void stats_construct(void) {
  int_stats.count = 0;
}

__attribute__ (( __destructor__ (101) ))
void stats_destruct(void) {
  size_t i;
  for (i = 0;i < int_stats.count; ++ i) {
    free(int_stats.names[i]);
  }
}

static void stats_print(FILE* out) {
  size_t i;
  for (i = 0; i < int_stats.count; ++ i) {
    fprintf(out, "%s = %d\n", int_stats.names[i], int_stats.values[i]);
  }
}

// Implementation of statistics
const statistics_ops_t statistics_ops = {
    stats_register_int,
    stats_print
};

