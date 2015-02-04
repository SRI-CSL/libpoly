/*
 * value.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

#include <value.h>

void value_construct(lp_value_t* v, lp_value_type_t type, const void* data);

void value_construct_copy(lp_value_t* v, const lp_value_t* from);

void value_destruct(lp_value_t* v);

void value_approx(const lp_value_t* v, interval_t* approx);

int value_print(const lp_value_t* v, FILE* out);

int value_cmp(const lp_value_t* v1, const lp_value_t* v2);

int value_cmp_void(const void* v1, const void* v2);
