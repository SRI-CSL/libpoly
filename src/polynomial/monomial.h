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

/** Construct an empty monomial */
void monomial_construct(const polynomial_context_t* ctx, monomial_t* m);

/** Construct an ordered copy of the monomial */
void monomial_construct_copy(const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort);

/** Destruct the monomial */
void monomial_destruct(monomial_t* m);

/** Clear the monomial to 0 */
void monomial_clear(const polynomial_context_t* ctx, monomial_t* m);

/** Assign another monomial */
void monomial_assign(const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort);

/** Add a variable power to the end of the monomial */
void monomial_push(monomial_t* m, variable_t x, unsigned d);

/** Remove a variable powerfrom the end of the monomial */
void monomial_pop(monomial_t* m);

/** Get the gcd of two monomials */
void monomial_gcd(const polynomial_context_t* ctx, monomial_t* gcd, const monomial_t* m1, const monomial_t* m2);

