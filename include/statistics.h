/*
 * statistics.h
 *
 *  Created on: Dec 4, 2013
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>

/**
 * Operations for management of statistics.
 */
typedef struct {

  /** Register an integer statistic */
  int* (*register_int) (const char* name);

  /** Output the statistics to a file */
  void (*print) (FILE* out);

} statistics_ops_t;

const statistics_ops_t statistics_ops;
