/*
 * monomial.h
 *
 *  Created on: Feb 10, 2014
 *      Author: dejan
 */

#pragma once

#include "number/integer.h"

#include "variable.h"

#include <stdio.h>

/**
 * Power of single variable.
 */
typedef struct {
  variable_t x;
  unsigned d;
} power_t;

/** Context for the polynomial operations */
typedef struct polynomial_context_struct polynomial_context_t;

/**
 * A monomial of the form a*x_1^d_1*...*x_size*d_size.
 */
typedef struct {
  /** The coefficient */
  integer_t a;
  /** Number of variables */
  size_t n;
  /** The capacity of the power array */
  size_t capacity;
  /** Array of variable powers */
  power_t* p;
} monomial_t;

typedef struct {
  /** Construct an empty monomial */
  void (*construct) (const polynomial_context_t* ctx, monomial_t* m);
  /** Construct an ordered copy of the monomial */
  void (*construct_copy) (const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort);
  /** Destruct the monomial */
  void (*destruct) (monomial_t* m);
  /** Clear the monomial to 0 */
  void (*clear) (const polynomial_context_t* ctx, monomial_t* m);
  /** Assign another monomial */
  void (*assign) (const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort);
  /** Print the monomial */
  int (*print) (const polynomial_context_t* ctx, const monomial_t* m, FILE* out);
  /** Add a variable power to the end of the monomial */
  void (*push) (monomial_t* m, variable_t x, unsigned d);
  /** Remove a variable powerfrom the end of the monomial */
  void (*pop) (monomial_t* m);
  /** Get the gcd of two monomials */
  void (*gcd) (const polynomial_context_t* ctx, monomial_t* gcd, const monomial_t* m1, const monomial_t* m2);

} monomial_ops_t;

/** Implementation of the monomial operations */
extern const monomial_ops_t monomial_ops;
