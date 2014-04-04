/*
 * output.h
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#pragma once

#include <output_language.h>

extern poly_output_language_t output_language;

void set_output_language(poly_output_language_t lang);

const char* get_power_symbol(void);
