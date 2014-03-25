/*
 * interval_arithmetic.h
 *
 *  Created on: Mar 24, 2014
 *      Author: dejan
 */

#include "interval/interval.h"

void interval_add(interval_t* S, const interval_t* I1, const interval_t* I2);

void interval_neg(interval_t* N, const interval_t* I);

void interval_sub(interval_t* S, const interval_t* I1, const interval_t* I2);

void interval_mul(interval_t* P, const interval_t* I1, const interval_t* I2);

void interval_pow(interval_t* P, const interval_t* I, unsigned n);

void dyadic_interval_add(dyadic_interval_t* S, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

void dyadic_interval_neg(dyadic_interval_t* N, const dyadic_interval_t* I);

void dyadic_interval_sub(dyadic_interval_t* S, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

void dyadic_interval_mul(dyadic_interval_t* P, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

void dyadic_interval_pow(dyadic_interval_t* P, const dyadic_interval_t* I, unsigned n);
