/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <rational_interval.h>
#include <dyadic_interval.h>

void rational_interval_add(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2);

void rational_interval_neg(lp_rational_interval_t* N, const lp_rational_interval_t* I);

void rational_interval_sub(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2);

void rational_interval_mul(lp_rational_interval_t* P, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2);

void rational_interval_pow(lp_rational_interval_t* P, const lp_rational_interval_t* I, unsigned n);

void dyadic_interval_add(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_neg(lp_dyadic_interval_t* N, const lp_dyadic_interval_t* I);

void dyadic_interval_sub(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_mul(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

void dyadic_interval_pow(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I, unsigned n);

/** Compute an overapproximation of the n th positive root of the reals in
    I. The precision of the approximation is at least prec (roughly the
    overapproximation of the result is at least 2^(-prec/n) ) .
 */
void dyadic_interval_root_overapprox(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I, unsigned n, unsigned prec);
