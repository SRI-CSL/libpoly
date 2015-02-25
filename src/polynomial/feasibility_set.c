/*
 * feasibility_set.c
 *
 *  Created on: Feb 24, 2015
 *      Author: dejan
 */

#include <feasibility_set.h>
#include <stdlib.h>

struct feasibility_set_struct {

};

lp_feasibility_set_t* lp_feasibility_set_new() {
  return 0;
}

void lp_feasibility_set_delete(lp_feasibility_set_t* set) {
  free(set);
}

int lp_feasibility_set_is_empty(const lp_feasibility_set_t* set) {
  return set == 0;
}

lp_feasibility_set_t* lp_polynomial_get_feasible_set(const lp_polynomial_t* p, lp_sign_condition_t sgn_condition, const lp_feasibility_set_t* domain) {
  (void)p;
  (void)sgn_condition;
  (void)domain;
  return 0;
}

int lp_feasibility_set_contains(const lp_feasibility_set_t* set, const lp_value_t* value) {
  (void)set;
  (void)value;
  return 1;
}
