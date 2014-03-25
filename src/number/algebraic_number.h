/*
 * algebraic_number.h
 *
 *  Created on: Mar 25, 2014
 *      Author: dejan
 */

#pragma once

#include <algebraic_number.h>

/**
 * Construct the algebraic number given it's polynomial and the isolating
 * interval. The number takes over the reference of f.
 */
void algebraic_number_construct(algebraic_number_t* a, upolynomial_t* f, const dyadic_interval_t* I);

/** Construct a zero algebraic number */
void algebraic_number_construct_zero(algebraic_number_t* a);

/** Construct a copy of the algebraic number. */
void algebraic_number_construct_copy(algebraic_number_t* a1, const algebraic_number_t* a2);

/** Construct the algebraic number from a dyadic rational */
void algebraic_number_construct_from_dyadic_rational(algebraic_number_t* a, const dyadic_rational_t* q);

/** Destruct the number */
void algebraic_number_destruct(algebraic_number_t* a);

/** Compare two algebraic numbers */
int algebraic_number_cmp(const algebraic_number_t* a1, const algebraic_number_t* a2);

/** Void version of the comparison, use with care. */
int algebraic_number_cmp_void(const void* a1, const void* a2);

/** Print the number */
int algebraic_number_print(const algebraic_number_t* a, FILE* out);

/** Return a string representation of the number */
char* algebraic_number_to_string(const algebraic_number_t* a);

/** Convert to double */
double algebraic_number_to_double(const algebraic_number_t* a);

/** Refine the number by halfing it's interval. */
void algebraic_number_refine(algebraic_number_t* a);

/** Addition */
void algebraic_number_add(algebraic_number_t* sum, const algebraic_number_t* a, const algebraic_number_t* b);

/** Subtraction */
void algebraic_number_sub(algebraic_number_t* sub, const algebraic_number_t* a, const algebraic_number_t* b);

/** Multiplication */
void algebraic_number_mul(algebraic_number_t* mul, const algebraic_number_t* a, const algebraic_number_t* b);

/** Multiplication */
void algebraic_number_pow(algebraic_number_t* pow, const algebraic_number_t* a, unsigned n);
