/*
 * feasibility_set.h
 *
 *  Created on: Feb 24, 2015
 *      Author: dejan
 */

#pragma once

#include "poly.h"
#include "sign_condition.h"

/**
 * Create a new feasibility set [-inf, inf].
 */
lp_feasibility_set_t* lp_feasibility_set_new();

/**
 * Delete the given feasibility set.
 */
void lp_feasibility_set_delete(lp_feasibility_set_t* set);

/**
 * Check if the given set is empty.
 */
int lp_feasibility_set_is_empty(const lp_feasibility_set_t* set);

/**
 * Check if the given value belongs to the set.
 */
int lp_feasibility_set_contains(const lp_feasibility_set_t* set, const lp_value_t* value);

/**
 * Pick a value from the feasible set (must be non-empty).
 */
lp_value_t* lp_feasibility_set_pick_value(const lp_feasibility_set_t* set);

/**
 * Print the set.
 */
int lp_feasibility_set_print(const lp_feasibility_set_t* set, FILE* out);

