/*
 * gcd.h
 *
 *  Created on: Mar 28, 2014
 *      Author: dejan
 */

#pragma once

#include "polynomial/coefficient.h"

/** Compute the greatest common divisor gcd(C1, C2). */
void coefficient_gcd(const lp_polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the primitive part pp(A(divide with content).
 */
void coefficient_pp(const lp_polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C);

/**
 * Compute the content cont(A(gcd of coefficients). Content is computed so that
 * A/cont(A) has a positive leading coefficient.
 */
void coefficient_cont(const lp_polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C);

/** Comput pp and cont at the same time */
void coefficient_pp_cont(const lp_polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C);

/** Compute the least common multiple. */
void coefficient_lcm(const lp_polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2);


