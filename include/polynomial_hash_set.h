/*
 * lp_polynomial_hash_set.h
 *
 *  Created on: May 27, 2015
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>

struct lp_polynomial_hash_set_struct {
  /** The data */
  lp_polynomial_t** data;
  /** Size of the data */
  size_t data_size;
  /** Number of set elements */
  size_t size;
  /** Treshold for resize */
  size_t resize_threshold;
};

/** Construct a new set */
void lp_polynomial_hash_set_construct(lp_polynomial_hash_set_t* set);

/** Destruct the set */
void lp_polynomial_hash_set_destruct(lp_polynomial_hash_set_t* set);

/** Returns true if empty */
int lp_polynomial_hash_set_is_empty(lp_polynomial_hash_set_t* set);

/** Check whether p is in set. */
int lp_polynomial_hash_set_contains(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p);

/** Add polynomial p to set. Returns true if p was added (not already in the set). */
int lp_polynomial_hash_set_insert(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p);

/** Close the set: compact the data so that all elements get stored in data[0..size]. No addition after close! */
void lp_polynomial_hash_set_close(lp_polynomial_hash_set_t* set);

/** Clear the set. */
void lp_polynomial_hash_set_clear(lp_polynomial_hash_set_t* set);

/** Print the set. */
int lp_polynomial_hash_set_print(const lp_polynomial_hash_set_t* set, FILE* out);
