/*
 * rational.c
 *
 *Created on: Jan 13, 2014
 *Author: dejan
 */

#include "rational.h"
#include "rational_internal.h"

void rational_register_prinf_extension(void) {
  assert(0);
}

const rational_ops_t rational_ops = {
    rational_construct,
    rational_construct_from_int,
    rational_construct_from_integer,
    rational_construct_from_double,
    rational_construct_from_dyadic,
    rational_construct_copy,
    rational_assign,
    rational_assign_int,
    rational_destruct,
    rational_print,
    rational_to_string,
    rational_to_double,
    rational_sgn,
    rational_cmp,
    rational_swap,
    rational_add,
    rational_add_integer,
    rational_sub,
    rational_neg,
    rational_inv,
    rational_mul,
    rational_mul_2exp,
    rational_pow,
    rational_div,
    rational_div_2exp,
    rational_register_prinf_extension,
};
