/*
 * output.h
 *
 *  Created on: Apr 4, 2014
 *      Author: dejan
 */

#include "upolynomial/umonomial.h"
#include "upolynomial/upolynomial.h"
#include "upolynomial/upolynomial_dense.h"
#include "upolynomial/factors.h"

#include <stdio.h>

int umonomial_print(const umonomial_t* m, FILE* out);

int upolynomial_print(const upolynomial_t* p, FILE* out);

char* upolynomial_to_string(const upolynomial_t* p);

int upolynomial_factors_print(const upolynomial_factors_t* f, FILE* out);

int upolynomial_dense_print(const upolynomial_dense_t* p_d, FILE* file);
