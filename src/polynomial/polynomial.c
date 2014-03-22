/*
 * polynomial.c
 *
 *  Created on: Jan 30, 2014
 *      Author: dejan
 */

#include "polynomial/polynomial_internal.h"

#include <printf.h>

static
int polynomial_printf_handler_set = 0;

static int
polynomial_printf_handler(FILE* out, const struct printf_info* info, const void *const *args) {
  const polynomial_t* P;
  // Get the coefficient
  P = *((const polynomial_t **) (args[0]));
  // Print the polynomial
  return polynomial_ops.print(P, out);
}

static int
polynomial_printf_arginfo(const struct printf_info *info, size_t n, int *argtypes, int *size) {
  if (n > 0) {
    argtypes[0] = PA_POINTER;
  }
  return 1;
}

static
void polynomial_register_printf_extension(void) {
  if (!polynomial_printf_handler_set) {
    polynomial_printf_handler_set = 1;
    register_printf_specifier('P', polynomial_printf_handler, polynomial_printf_arginfo);
  }
}

/** Set the power symbol for printouts */
void polynomial_set_power_symbol(const char* power) {
  coefficient_ops.set_power_symbol(power);
}


const polynomial_ops_t polynomial_ops = {
  polynomial_construct,
  polynomial_construct_simple,
  polynomial_construct_copy,
  polynomial_destruct,
  polynomial_alloc,
  polynomial_new,
  polynomial_swap,
  polynomial_assign,
  polynomial_context,
  polynomial_degree,
  polynomial_top_variable,
  polynomial_get_coefficient,
  polynomial_is_constant,
  polynomial_is_zero,
  polynomial_sgn,
  polynomial_cmp,
  polynomial_cmp_type,
  polynomial_divides,
  polynomial_print,
  polynomial_to_string,
  polynomial_add,
  polynomial_sub,
  polynomial_neg,
  polynomial_mul,
  polynomial_shl,
  polynomial_pow,
  polynomial_add_mul,
  polynomial_sub_mul,
  polynomial_reduce,
  polynomial_div,
  polynomial_rem,
  polynomial_divrem,
  polynomial_derivative,
  polynomial_gcd,
  polynomial_lcm,
  polynomial_resultant,
  polynomial_psc,
  polynomial_set_power_symbol,
  polynomial_register_printf_extension
};

