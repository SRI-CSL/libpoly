/*
 * assignment_internal.h
 *
 *  Created on: Mar 7, 2014
 *      Author: dejan
 */

#pragma once

#include <assignment.h>

void value_construct(value_t* v, value_type_t type, const void* data);

void value_construct_copy(value_t* v, const value_t* from);

void value_destruct(value_t* v);

void value_approx(const value_t* v, interval_t* approx);

int value_print(const value_t* v, FILE* out);

void assignment_construct(assignment_t* m, const variable_db_t* var_db);

assignment_t* assignment_new(const variable_db_t* var_db);

void assignment_destruct(assignment_t* m);

void assignment_delete(assignment_t* m);

int assignment_print(const assignment_t* m, FILE* out);

char* assignment_to_string(const assignment_t* m);

void assignment_set_value(assignment_t* m, variable_t x, const value_t* value);

const value_t* assignment_get_value(const assignment_t* m, variable_t x);

void assignment_get_value_approx(const assignment_t* m, variable_t x, interval_t* approx);

int assignment_sgn(const assignment_t* m, const polynomial_t* A);
