/*
 * value.c
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#include "number/value.h"

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"
#include "number/algebraic_number.h"
#include "interval/interval.h"

void value_construct(value_t* v, value_type_t type, const void* data) {
  v->type = type;
  switch(type) {
  case VALUE_NONE:
    break;
  case VALUE_INTEGER:
    integer_construct_copy(Z, &v->value.z, data);
    break;
  case VALUE_RATIONAL:
    rational_construct_copy(&v->value.q, data);
    break;
  case VALUE_DYADIC_RATIONAL:
    dyadic_rational_construct_copy(&v->value.dy_q, data);
    break;
  case VALUE_ALGEBRAIC:
    algebraic_number_construct_copy(&v->value.a, data);
    break;
  }
}

void value_construct_copy(value_t* v, const value_t* from) {
  switch(from->type) {
  case VALUE_NONE:
    value_construct(v, VALUE_NONE, 0);
    break;
  case VALUE_INTEGER:
    value_construct(v, VALUE_INTEGER, &from->value.z);
    break;
  case VALUE_RATIONAL:
    value_construct(v, VALUE_RATIONAL, &from->value.q);
    break;
  case VALUE_DYADIC_RATIONAL:
    value_construct(v, VALUE_DYADIC_RATIONAL, &from->value.dy_q);
    break;
  case VALUE_ALGEBRAIC:
    value_construct(v, VALUE_ALGEBRAIC, &from->value.a);
    break;
  }
}

void value_destruct(value_t* v) {
  switch(v->type) {
  case VALUE_NONE:
    break;
  case VALUE_INTEGER:
    integer_destruct(&v->value.z);
    break;
  case VALUE_RATIONAL:
    rational_destruct(&v->value.q);
    break;
  case VALUE_DYADIC_RATIONAL:
    dyadic_rational_destruct(&v->value.dy_q);
    break;
  case VALUE_ALGEBRAIC:
    algebraic_number_destruct(&v->value.a);
    break;
  }
}

void value_approx(const value_t* v, interval_t* approx) {
  switch (v->type) {
  case VALUE_INTEGER:
    assert(0);
    break;
  case VALUE_RATIONAL:
    interval_construct(approx, &v->value.q, 0, &v->value.q, 0);
    break;
  case VALUE_DYADIC_RATIONAL:
    interval_construct_from_dyadic(approx, &v->value.dy_q, 0, &v->value.dy_q, 0);
    break;
  case VALUE_ALGEBRAIC:
    assert(0);
    break;
  case VALUE_NONE:
    assert(0);
    break;
  }
}

int value_print(const value_t* v, FILE* out) {
  int ret = 0;
  switch (v->type) {
  case VALUE_NONE:
    ret += fprintf(out, "<null>");
    break;
  case VALUE_INTEGER:
    ret += integer_print(&v->value.z, out);
    break;
  case VALUE_RATIONAL:
    ret += rational_print(&v->value.q, out);
    break;
  case VALUE_DYADIC_RATIONAL:
    ret += dyadic_rational_print(&v->value.dy_q, out);
    break;
  case VALUE_ALGEBRAIC:
    ret += algebraic_number_print(&v->value.a, out);
    break;
  }
  return ret;
}

int value_cmp(const value_t* v1, const value_t* v2) {
  return v1 == v2;
}

int value_cmp_void(const void* v1, const void* v2) {
  return v1 == v2;
}

const value_ops_t value_ops = {
    value_construct,
    value_construct_copy,
    value_destruct,
    value_approx,
    value_cmp,
    value_cmp_void,
    value_print
};
