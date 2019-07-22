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

#include "statistics.h"

#include <string.h>
#include <stdlib.h>

#define INITIAL_STATS 100

typedef struct int_stats_struct {
  size_t count;
  int values[INITIAL_STATS];
  char* names[INITIAL_STATS];
} int_stats_t;

/** Integer statistics (default init all to 0) */
__thread
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

__attribute__ (( __destructor__ ))
void stats_destruct(void) {
  size_t i;
  for (i = 0;i < int_stats.count; ++ i) {
    free(int_stats.names[i]);
  }
}
