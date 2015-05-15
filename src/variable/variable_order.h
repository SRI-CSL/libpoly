/*
 * variable_order.h
 *
 *  Created on: May 15, 2015
 *      Author: dejan
 */

#pragma once

#include <poly.h>

/** Make a variable the top variable (bigger than anything) */
void lp_variable_order_make_top(lp_variable_order_t* var_order, lp_variable_t var);

/** Make a variable the bottom variable (smaller than anything) */
void lp_variable_rder_make_bot(lp_variable_order_t* var_order, lp_variable_t var);
