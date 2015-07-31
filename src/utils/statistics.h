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

/** Register a new statistic with the given name */
int* stats_register_int(const char* name);

/** Print the statistics to the given file */
void stats_print(FILE* out);

// names and output
#define STAT_NAME(module, name) __stat_ ## module ## _ ## name
#define STAT_INIT_NAME(module, name) __stat_ ## module ## _ ## name ## __init
#define STAT_OUT(module, name) #module "::" #name

/**
 * Use to declare a statistic of the given type with the variable name
 * c_name and the output string string_name. The statistic will be registered
 * at library load_time. This should only be used at global scope in the .c
 * files.
 */
#define STAT_DECLARE(type, module, name)      \
  type* STAT_NAME(module, name);              \
                                              \
  __attribute__ (( __constructor__ ))   \
  void STAT_INIT_NAME(module, name) (void) {  \
    STAT_NAME(module, name) = stats_register_ ## type(STAT_OUT(module, name)); \
  }

/**
 * Use to reference the (previously declared) statistic.
 */
#define STAT(module, name) (*STAT_NAME(module, name))
