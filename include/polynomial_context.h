/*
 * polynomial_context.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "variable.h"
#include "integer.h"
#include "variable_order.h"

/**
 * Information providing semantic context for available operations on
 * polynomials.
 */
struct polynomial_context_struct {
  /** Reference count of this polynomial */
  size_t ref_count;
  /** The ring of base operations */
  int_ring_t* K;
  /** The variable database */
  variable_db_t* var_db;
  /** The order of variables */
  variable_order_t* var_order;
};

/** Operations on the polynomial context */
typedef struct {

  /** Create a new context and attach. */
  polynomial_context_t* (*new) (int_ring_t* K, variable_db_t* var_db, variable_order_t* var_order);

  /** Attach to existing context */
  void (*attach) (polynomial_context_t* ctx);

  /** Delete the context (dec reference count and if 0 also free the memory) */
  void (*detach) (polynomial_context_t* ctx);

  /**
   * Compare two contexts. Two contexts are equal if their rings test for
   * equality and their variable order is identical (pointers)
   */
  int (*equal) (const polynomial_context_t* ctx1, const polynomial_context_t* ctx2);

} polynomial_context_ops_t;

/** Implementation of the context operations */
extern const polynomial_context_ops_t polynomial_context_ops;
