/*
 * output.h
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#include <upolynomial.h>
#include <upolynomial_factors.h>

#include "upolynomial/umonomial.h"
#include "upolynomial/upolynomial_dense.h"

#include <stdio.h>

int umonomial_print(const ulp_monomial_t* m, FILE* out);

int upolynomial_dense_print(const upolynomial_dense_t* p_d, FILE* file);
