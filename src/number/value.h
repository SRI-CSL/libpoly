/*
 * value.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

#include <value.h>

void lp_value_construct(lp_value_t* v, lp_value_type_t type, const void* data);

void lp_value_construct_copy(lp_value_t* v, const lp_value_t* from);

void lp_value_destruct(lp_value_t* v);

void lp_value_approx(const lp_value_t* v, lp_interval_t* approx);

int lp_value_print(const lp_value_t* v, FILE* out);

int lp_value_cmp(const lp_value_t* v1, const lp_value_t* v2);

int lp_value_cmp_void(const void* v1, const void* v2);
