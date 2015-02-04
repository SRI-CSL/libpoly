/*
 * assignment_internal.h
 *
 *  Created on: Mar 7, 2014
 *      Author: dejan
 */

#pragma once

#include <assignment.h>

void value_construct(lp_value_t* v, lp_value_type_t type, const void* data);

void value_construct_copy(lp_value_t* v, const lp_value_t* from);

void value_destruct(lp_value_t* v);

void value_approx(const lp_value_t* v, interval_t* approx);

int value_print(const lp_value_t* v, FILE* out);

void assignment_construct(lp_assignment_t* m, const lp_variable_db_t* var_db);

lp_assignment_t* assignment_new(const lp_variable_db_t* var_db);

void assignment_destruct(lp_assignment_t* m);

void assignment_delete(lp_assignment_t* m);

int assignment_print(const lp_assignment_t* m, FILE* out);

char* assignment_to_string(const lp_assignment_t* m);

void assignment_set_value(lp_assignment_t* m, lp_variable_t x, const lp_value_t* value);

const lp_value_t* assignment_get_value(const lp_assignment_t* m, lp_variable_t x);

void assignment_get_value_approx(const lp_assignment_t* m, lp_variable_t x, interval_t* approx);

int assignment_sgn(const lp_assignment_t* m, const lp_polynomial_t* A);
