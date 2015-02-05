/*
 * polynomial_context.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dejan
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
  /** Reference count of this polynomial */
  size_t ref_count;
  /** The ring of base operations */
  lp_int_ring_t* K;
  /** The variable database */
  lp_variable_db_t* var_db;
  /** The order of variables */
  lp_variable_order_t* var_order;
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
