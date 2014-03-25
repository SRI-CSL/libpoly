/*
 * interval_internal.h
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#pragma once

#include <interval.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

#include <assert.h>
#include <limits.h>

// Rational interval

void interval_construct(interval_t* I, const rational_t* a, int a_open, const rational_t* b, int b_open);

void interval_construct_point(interval_t* I, const rational_t* a);

void interval_construct_zero(interval_t* I);

void interval_construct_copy(interval_t* I, const interval_t* from);

void interval_construct_from_dyadic(interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open);

void interval_construct_from_int(interval_t* I, long a, int a_open, long b, int b_open);

void interval_construct_from_integer(interval_t* I, const integer_t* a, int a_open, const integer_t* b, int b_open);

void interval_construct_point(interval_t* I, const rational_t* q);

void interval_destruct(interval_t* I);

void interval_assign(interval_t* I, const interval_t* from);

void interval_swap(interval_t* I1, interval_t* I2);

int interval_sgn(const interval_t* I);

int interval_print(const interval_t* I, FILE* out);

// Dyadic rational

void dyadic_interval_construct(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open);

void dyadic_interval_construct_zero(dyadic_interval_t* I);

void dyadic_interval_construct_copy(dyadic_interval_t* I, const dyadic_interval_t* from);

void dyadic_interval_construct_from_dyadic(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open);

void dyadic_interval_construct_from_int(dyadic_interval_t* I, long a, int a_open, long b, int b_open);

void dyadic_interval_construct_from_integer(dyadic_interval_t* I, const integer_t* a, int a_open, const integer_t* b, int b_open);

void dyadic_interval_construct_point(dyadic_interval_t* I, const dyadic_rational_t* q);

void dyadic_interval_construct_from_split(dyadic_interval_t* I_left, dyadic_interval_t* I_right, const dyadic_interval_t* I, int left_open, int right_open);

void dyadic_interval_destruct(dyadic_interval_t* I);

void dyadic_interval_assign(dyadic_interval_t* I, const dyadic_interval_t* from);

void dyadic_interval_swap(dyadic_interval_t* I1, dyadic_interval_t* I2);

int dyadic_interval_sgn(const dyadic_interval_t* I);

int dyadic_interval_contains(const dyadic_interval_t* I, const dyadic_rational_t* q);

void dyadic_interval_construct_intersection(dyadic_interval_t* I, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

void dyadic_interval_collapse_to(dyadic_interval_t* I, const dyadic_rational_t* q);

void dyadic_interval_set_a(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open);

void dyadic_interval_set_b(dyadic_interval_t* I, const dyadic_rational_t* b, int b_open);

int dyadic_interval_equals(const dyadic_interval_t* I1, const dyadic_interval_t* I2);

int dyadic_interval_disjunct(const dyadic_interval_t* I1, const dyadic_interval_t* I2);

void dyadic_interval_scale(dyadic_interval_t* I, int n);

int dyadic_interval_print(const dyadic_interval_t* I, FILE* out);

