/*
 * value.c
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#include "number/value.h"

#include <interval.h>
#include <algebraic_number.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

void lp_value_construct(lp_value_t* v, lp_value_type_t type, const void* data) {
  v->type = type;
  switch(type) {
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
  case LP_VALUE_NONE:
    break;
  case LP_VALUE_INTEGER:
    integer_construct_copy(lp_Z, &v->value.z, data);
    break;
  case LP_VALUE_RATIONAL:
    rational_construct_copy(&v->value.q, data);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_construct_copy(&v->value.dy_q, data);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_copy(&v->value.a, data);
    break;
  }
}

lp_value_t* lp_value_new(lp_value_type_t type, const void* data) {
  lp_value_t* result = malloc(sizeof(lp_value_t));
  lp_value_construct(result, type, data);
  return result;
}

void lp_value_construct_copy(lp_value_t* v, const lp_value_t* from) {
  switch(from->type) {
  case LP_VALUE_NONE:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    lp_value_construct(v, from->type, 0);
    break;
  case LP_VALUE_INTEGER:
    lp_value_construct(v, LP_VALUE_INTEGER, &from->value.z);
    break;
  case LP_VALUE_RATIONAL:
    lp_value_construct(v, LP_VALUE_RATIONAL, &from->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_value_construct(v, LP_VALUE_DYADIC_RATIONAL, &from->value.dy_q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_value_construct(v, LP_VALUE_ALGEBRAIC, &from->value.a);
    break;
  }
}

void lp_value_destruct(lp_value_t* v) {
  switch(v->type) {
  case LP_VALUE_NONE:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    break;
  case LP_VALUE_INTEGER:
    integer_destruct(&v->value.z);
    break;
  case LP_VALUE_RATIONAL:
    rational_destruct(&v->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_destruct(&v->value.dy_q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_destruct(&v->value.a);
    break;
  }
}

void lp_value_delete(lp_value_t* v) {
  lp_value_destruct(v);
  free(v);
}


void lp_value_approx(const lp_value_t* v, lp_interval_t* approx) {
  switch (v->type) {
  case LP_VALUE_INTEGER:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    assert(0);
    break;
  case LP_VALUE_RATIONAL:
    lp_interval_construct_point(approx, &v->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_interval_construct_from_dyadic(approx, &v->value.dy_q, 0, &v->value.dy_q, 0);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_interval_construct_from_dyadic_interval(approx, &v->value.a.I);
    assert(0);
    break;
  case LP_VALUE_NONE:
    assert(0);
    lp_interval_construct_zero(approx);
    break;
  }
}

int lp_value_print(const lp_value_t* v, FILE* out) {
  int ret = 0;
  switch (v->type) {
  case LP_VALUE_NONE:
    ret += fprintf(out, "<null>");
    break;
  case LP_VALUE_PLUS_INFINITY:
    ret += fprintf(out, "+inf");
    break;
  case LP_VALUE_MINUS_INFINITY:
    ret += fprintf(out, "-inf");
    break;
  case LP_VALUE_INTEGER:
    ret += integer_print(&v->value.z, out);
    break;
  case LP_VALUE_RATIONAL:
    ret += rational_print(&v->value.q, out);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    ret += dyadic_rational_print(&v->value.dy_q, out);
    break;
  case LP_VALUE_ALGEBRAIC:
    ret += lp_algebraic_number_print(&v->value.a, out);
    break;
  }
  return ret;
}

int lp_value_cmp(const lp_value_t* v1, const lp_value_t* v2) {
  assert(0);
  return v1 == v2;
}

int lp_value_cmp_void(const void* v1, const void* v2) {
  assert(0);
  return v1 == v2;
}
