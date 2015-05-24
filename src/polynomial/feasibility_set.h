/*
 * feasibility_set.h
 *
 *  Created on: May 21, 2015
 *      Author: dejan
 */

#pragma once

void lp_feasibility_set_construct(lp_feasibility_set_t* s, size_t size);

lp_feasibility_set_t* lp_feasibility_set_new_internal(size_t size);
