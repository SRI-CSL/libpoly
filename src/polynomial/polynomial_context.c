/*
 * polynomial_context.c
 *
 *  Created on: Feb 13, 2014
 *      Author: dejan
 */

#include <variable_db.h>
#include <polynomial_context.h>

#include <stdlib.h>
#include <assert.h>

#define TEMP_VARIABLE_SIZE 10

static
void lp_polynomial_context_construct(lp_polynomial_context_t* ctx, lp_int_ring_t* K, lp_variable_db_t* var_db, lp_variable_order_t* var_order) {
  ctx->ref_count = 0;
  ctx->var_db = var_db;
  ctx->K = K;
  ctx->var_db = var_db;
  ctx->var_order = var_order;

  ctx->var_tmp = malloc(sizeof(lp_variable_t)*TEMP_VARIABLE_SIZE);
  size_t i;
  for (i = 0; i < TEMP_VARIABLE_SIZE; ++ i) {
    char name[10];
    sprintf(name, "#%zu", i);
    ctx->var_tmp[i] = lp_variable_db_new_variable(var_db, name);
  }
}

static
void lp_polynomial_context_destruct(lp_polynomial_context_t* ctx) {
  free(ctx->var_tmp);
}

void lp_polynomial_context_attach(lp_polynomial_context_t* ctx) {
  if (ctx->K) {
    lp_int_ring_attach(ctx->K);
  }
  if (ctx->var_db) {
    lp_variable_db_attach(ctx->var_db);
  }
  if (ctx->var_order) {
    lp_variable_order_attach(ctx->var_order);
  }
  ctx->ref_count ++;
}

void lp_polynomial_context_detach(lp_polynomial_context_t* ctx) {
  if (ctx->K) {
    lp_int_ring_detach(ctx->K);
  }
  if (ctx->var_db) {
    lp_variable_db_detach(ctx->var_db);
  }
  if (ctx->var_order) {
    lp_variable_order_detach(ctx->var_order);
  }
  assert(ctx->ref_count > 0);
  ctx->ref_count --;
  if (ctx->ref_count == 0) {
    lp_polynomial_context_destruct(ctx);
  }
}

lp_polynomial_context_t* lp_polynomial_context_new(lp_int_ring_t* K, lp_variable_db_t* var_db, lp_variable_order_t* var_order) {
  lp_polynomial_context_t* result = malloc(sizeof(lp_polynomial_context_t));
  lp_polynomial_context_construct(result, K, var_db, var_order);
  lp_polynomial_context_attach(result);
  return result;
}


int lp_polynomial_context_equal(const lp_polynomial_context_t* ctx1, const lp_polynomial_context_t* ctx2) {
  if (ctx1 == ctx2) return 1;
  if (ctx1 && ctx2) {
    return lp_int_ring_equal(ctx1->K, ctx2->K) && ctx1->var_order == ctx2->var_order;
  } else {
    return 0;
  }
}

/** Get a temp variable */
lp_variable_t lp_polynomial_context_get_temp_variable(const lp_polynomial_context_t* ctx_const) {
  lp_polynomial_context_t* ctx = (lp_polynomial_context_t*) ctx_const;
  assert(ctx->var_tmp_size < TEMP_VARIABLE_SIZE);
  return ctx->var_tmp[ctx->var_tmp_size ++];
}

/** Release the variable (has to be the last one obtained and not released */
void lp_polynomial_context_release_temp_variable(const lp_polynomial_context_t* ctx_const, lp_variable_t x) {
  lp_polynomial_context_t* ctx = (lp_polynomial_context_t*) ctx_const;
  assert(ctx->var_tmp_size > 0);
  ctx->var_tmp_size --;
  assert(ctx->var_tmp[ctx->var_tmp_size] == x);
}
