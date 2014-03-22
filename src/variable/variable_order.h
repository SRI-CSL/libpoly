/*
 * variable_order.h
 *
 *  Created on: Mar 22, 2014
 *      Author: dejan
 */

#pragma once

#include <variable_order.h>

void variable_order_simple_construct(variable_order_simple_t* var_order);

void variable_order_simple_destruct(variable_order_simple_t* var_order);

void variable_order_simple_attach(variable_order_t* var_order);

void variable_order_simple_detach(variable_order_t* var_order);

variable_order_t* variable_order_simple_new(void);

int variable_order_simple_cmp(const variable_order_t* var_order, variable_t x, variable_t y);

size_t variable_order_simple_size(const variable_order_simple_t* var_order);

void variable_order_simple_clear(variable_order_simple_t* var_order);

void variable_order_simple_push(variable_order_simple_t* var_order, variable_t var);

void variable_order_simple_pop(variable_order_simple_t* var_order);

int variable_order_simple_print(const variable_order_simple_t* var_order, const variable_db_t* var_db, FILE* out);

char* variable_order_simple_to_string(const variable_order_simple_t* var_order, const variable_db_t* var_db);

int variable_order_simple_contains(variable_order_simple_t* var_order, variable_t x);
