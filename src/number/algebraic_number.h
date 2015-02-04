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
void algebraic_number_construct(lp_algebraic_number_t* a, lp_upolynomial_t* f, const lp_dyadic_interval_t* I);

/** Construct a zero algebraic number */
void algebraic_number_construct_zero(lp_algebraic_number_t* a);

/** Construct a copy of the algebraic number. */
void algebraic_number_construct_copy(lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

/** Construct the algebraic number from a dyadic rational */
void algebraic_number_construct_from_dyadic_rational(lp_algebraic_number_t* a, const lp_dyadic_rational_t* q);

/** Destruct the number */
void algebraic_number_destruct(lp_algebraic_number_t* a);

/** Swap the two numbers */
void algebraic_number_swap(lp_algebraic_number_t* a, lp_algebraic_number_t* b);

/** Compare two algebraic numbers */
int algebraic_number_cmp(const lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

/** Void version of the comparison, use with care. */
int algebraic_number_cmp_void(const void* a1, const void* a2);

/** Print the number */
int algebraic_number_print(const lp_algebraic_number_t* a, FILE* out);

/** Return a string representation of the number */
char* algebraic_number_to_string(const lp_algebraic_number_t* a);

/** Convert to double */
double algebraic_number_to_double(const lp_algebraic_number_t* a);

/** Refine the number by halfing it's interval. */
void algebraic_number_refine(lp_algebraic_number_t* a);

/** Addition */
void algebraic_number_add(lp_algebraic_number_t* sum, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Subtraction */
void algebraic_number_sub(lp_algebraic_number_t* sub, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Negation */
void algebraic_number_neg(lp_algebraic_number_t* neg, const lp_algebraic_number_t* a);

/** Multiplication */
void algebraic_number_mul(lp_algebraic_number_t* mul, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Multiplication */
void algebraic_number_pow(lp_algebraic_number_t* pow, const lp_algebraic_number_t* a, unsigned n);
