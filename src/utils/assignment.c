/*
 * assignment.c
 *
 *  Created on: Feb 25, 2014
 *      Author: dejan
 */

#include "assignment.h"
#include "assignment_internal.h"

#include <malloc.h>

const value_ops_t value_ops = {
    value_construct,
    value_construct_copy,
    value_destruct,
    value_approx,
    value_print
};

const assignment_ops_t assignment_ops = {
    assignment_construct,
    assignment_new,
    assignment_destruct,
    assignment_print,
    assignment_to_string,
    assignment_set_value,
    assignment_get_value,
    assignment_get_value_approx,
    assignment_sgn
};
