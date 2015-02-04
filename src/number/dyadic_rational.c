/*
 * dyadic_dyadic_rational.c
 *
 *  Created on: Jan 22, 2014
 *      Author: dejan
 */

#include "number/dyadic_rational.h"

const lp_dyadic_rational_ops_t lp_dyadic_rational_ops = {
    dyadic_rational_construct,
    dyadic_rational_construct_from_int,
    dyadic_rational_construct_from_integer,
    dyadic_rational_construct_from_double,
    dyadic_rational_construct_copy,
    dyadic_rational_assign,
    dyadic_rational_assign_int,
    dyadic_rational_destruct,
    dyadic_rational_print,
    dyadic_rational_to_string,
    dyadic_rational_to_double,
    dyadic_rational_sgn,
    dyadic_rational_cmp,
    dyadic_rational_swap,
    dyadic_rational_add,
    dyadic_rational_add_integer,
    dyadic_rational_sub,
    dyadic_rational_neg,
    dyadic_rational_mul,
    dyadic_rational_mul_2exp,
    dyadic_rational_pow,
    dyadic_rational_div_2exp
};
