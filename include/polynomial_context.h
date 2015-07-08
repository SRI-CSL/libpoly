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

#include "poly.h"

#include "integer.h"
#include "variable_order.h"

/**
 * Information providing semantic context for available operations on
 * polynomials.
 */
struct lp_polynomial_context_struct {
  /** Reference count of this context */
  size_t ref_count;
  /** The ring of base operations */
  lp_int_ring_t* K;
  /** The variable database */
  lp_variable_db_t* var_db;
  /** The order of variables */
  lp_variable_order_t* var_order;
  /** Temporary variables for internal purposes */
  lp_variable_t* var_tmp;
  /** Size of temporary variables */
  size_t var_tmp_size;
};

/** Create a new context and attach. */
lp_polynomial_context_t* lp_polynomial_context_new(lp_int_ring_t* K, lp_variable_db_t* var_db, lp_variable_order_t* var_order);

/** Attach to existing context */
void lp_polynomial_context_attach(lp_polynomial_context_t* ctx);

/** Delete the context (dec reference count and if 0 also free the memory) */
void lp_polynomial_context_detach(lp_polynomial_context_t* ctx);

/**
 * Compare two contexts. Two contexts are equal if their rings test for
 * equality and their variable order is identical (pointers)
 */
int lp_polynomial_context_equal(const lp_polynomial_context_t* ctx1, const lp_polynomial_context_t* ctx2);
