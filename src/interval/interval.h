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

void lp_interval_construct(lp_interval_t* I, const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open);

void lp_interval_construct_point(lp_interval_t* I, const lp_rational_t* a);

void lp_interval_construct_zero(lp_interval_t* I);

void lp_interval_construct_copy(lp_interval_t* I, const lp_interval_t* from);

void lp_interval_construct_from_dyadic(lp_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

void lp_interval_construct_from_dyadic_interval(lp_interval_t* I, const lp_dyadic_interval_t* I_d);

void lp_interval_construct_from_int(lp_interval_t* I, long a, int a_open, long b, int b_open);

void lp_interval_construct_from_integer(lp_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

void lp_interval_construct_point(lp_interval_t* I, const lp_rational_t* q);

void lp_interval_destruct(lp_interval_t* I);

void lp_interval_assign(lp_interval_t* I, const lp_interval_t* from);

void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2);

int lp_interval_sgn(const lp_interval_t* I);

int lp_interval_contains_zero(const lp_interval_t* I);

int lp_interval_contains(const lp_interval_t* I, const lp_rational_t* q);

int lp_interval_print(const lp_interval_t* I, FILE* out);

// Dyadic rational

void lp_dyadic_interval_construct(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

void lp_dyadic_interval_construct_zero(lp_dyadic_interval_t* I);

void lp_dyadic_interval_construct_copy(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from);

void lp_dyadic_interval_construct_from_dyadic(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

void lp_dyadic_interval_construct_from_int(lp_dyadic_interval_t* I, long a, int a_open, long b, int b_open);

void lp_dyadic_interval_construct_from_integer(lp_dyadic_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

void lp_dyadic_interval_construct_point(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

void lp_dyadic_interval_construct_from_split(lp_dyadic_interval_t* I_left, lp_dyadic_interval_t* I_right, const lp_dyadic_interval_t* I, int left_open, int right_open);

void lp_dyadic_interval_destruct(lp_dyadic_interval_t* I);

void lp_dyadic_interval_assign(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from);

void lp_dyadic_interval_swap(lp_dyadic_interval_t* I1, lp_dyadic_interval_t* I2);

int lp_dyadic_interval_sgn(const lp_dyadic_interval_t* I);

int lp_dyadic_interval_contains(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

int lp_dyadic_interval_contains_zero(const lp_dyadic_interval_t* I);

void lp_dyadic_interval_construct_intersection(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void lp_dyadic_interval_collapse_to(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

void lp_dyadic_interval_set_a(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open);

void lp_dyadic_interval_set_b(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* b, int b_open);

int lp_dyadic_interval_equals(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

int lp_dyadic_interval_disjunct(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void lp_dyadic_interval_scale(lp_dyadic_interval_t* I, int n);

int lp_dyadic_interval_print(const lp_dyadic_interval_t* I, FILE* out);

