/*
 * value.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

#include <value.h>

void value_construct(value_t* v, value_type_t type, const void* data);

void value_construct_copy(value_t* v, const value_t* from);

void value_destruct(value_t* v);

void value_approx(const value_t* v, interval_t* approx);

int value_print(const value_t* v, FILE* out);

int value_cmp(const value_t* v1, const value_t* v2);

int value_cmp_void(const void* v1, const void* v2);
