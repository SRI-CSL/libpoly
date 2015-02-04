/*
 * interval_arithmetic.h
 *
 *  Created on: Mar 24, 2014
 *      Author: dejan
 */

#include <interval.h>

void interval_add(lp_interval_t* S, const lp_interval_t* I1, const lp_interval_t* I2);

void interval_neg(lp_interval_t* N, const lp_interval_t* I);

void interval_sub(lp_interval_t* S, const lp_interval_t* I1, const lp_interval_t* I2);

void interval_mul(lp_interval_t* P, const lp_interval_t* I1, const lp_interval_t* I2);

void interval_pow(lp_interval_t* P, const lp_interval_t* I, unsigned n);

void dyadic_interval_add(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_neg(lp_dyadic_interval_t* N, const lp_dyadic_interval_t* I);

void dyadic_interval_sub(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_mul(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_pow(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I, unsigned n);
