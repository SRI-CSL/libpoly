/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <polynomial_heap.h>
#include <polynomial.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

/** Default initial size (must be a power of 2) */
#define LP_POLYNOMIAL_HEAP_DEFAULT_SIZE 32

void lp_polynomial_heap_construct(lp_polynomial_heap_t *heap, lp_polynomial_heap_compare_f cmp) {
  heap->data = malloc(LP_POLYNOMIAL_HEAP_DEFAULT_SIZE * sizeof(lp_polynomial_t*));
  heap->data_size = LP_POLYNOMIAL_HEAP_DEFAULT_SIZE;
  heap->size = 0;
  heap->cmp = cmp;
}

void lp_polynomial_heap_destruct(lp_polynomial_heap_t* heap) {
  // remove all the polynomials
  for (size_t i = 0; i < heap->size; ++i) {
    lp_polynomial_delete(heap->data[i]);
  }
  // free the space
  free(heap->data);
  heap->data = NULL;
}

#define SWAP(A,B) ({ lp_polynomial_t *tmp = B; B = A; A = tmp; })
// macros assume a heap indices 1...size
#define HEAP_SWAP(heap, i, j) SWAP(heap->data[(i)-1], heap->data[(j)-1])
#define HEAP_CMP(heap, i, j) (heap->cmp(heap->data[(i)-1], heap->data[(j)-1]))

int lp_polynomial_heap_is_empty(lp_polynomial_heap_t *heap) {
  return heap->size == 0;
}

static
void lp_polynomial_heap_extend(lp_polynomial_heap_t *heap) {
  // double the size
  heap->data_size <<= 1;
  heap->data = realloc(heap->data, heap->data_size);
}

static
void lp_polynomial_heap_heapify_up(lp_polynomial_heap_t *heap) {
  for (size_t pos = heap->size;
       pos > 1 && HEAP_CMP(heap, pos / 2, pos) < 0;
       pos /= 2) {
    // if data[pos] is smaller or equal than data[parent] swap
    HEAP_SWAP(heap, pos / 2, pos);
  }
}

static
void lp_polynomial_heap_heapify_down(lp_polynomial_heap_t *heap) {
  size_t pos = 1;
  // while there is still at least a left child
  while (2 * pos <= heap->size) {
    size_t
        l = 2 * pos,
        r = 2 * pos + 1;

    // find larger child
    size_t o = (r <= heap->size && HEAP_CMP(heap, l, r) < 0) ? r : l;

    // in case current element is larger than larger child, we're done
    if (HEAP_CMP(heap, pos, o) >= 0) break;

    // swap larger child with current pos
    HEAP_SWAP(heap, pos, o);
    pos = o;
  }
}

void lp_polynomial_heap_push(lp_polynomial_heap_t* heap, const lp_polynomial_t* p) {
  heap->size++;
  if (heap->size > heap->data_size) {
    lp_polynomial_heap_extend(heap);
  }
  heap->data[heap->size - 1] = lp_polynomial_new_copy(p);
  lp_polynomial_heap_heapify_up(heap);
}

const lp_polynomial_t* lp_polynomial_heap_peek(lp_polynomial_heap_t* heap) {
  // empty heap does not have a top
  if (heap->size == 0) {
    return NULL;
  }
  // return the top element
  return heap->data[0];
}

lp_polynomial_t* lp_polynomial_heap_pop(lp_polynomial_heap_t* heap) {
  // empty heap does not have a top
  if (heap->size == 0) {
    return NULL;
  }
  // get the top element
  lp_polynomial_t *ret = heap->data[0];
  // reduce size and moves last element to top
  heap->data[0] = heap->data[--heap->size];
  // restore heap condition
  lp_polynomial_heap_heapify_down(heap);

  return ret;
}

void lp_polynomial_heap_clear(lp_polynomial_heap_t* heap) {
  lp_polynomial_heap_compare_f cmp = heap->cmp;
  lp_polynomial_heap_destruct(heap);
  lp_polynomial_heap_construct(heap, cmp);
}

void lp_polynomial_heap_print(lp_polynomial_heap_t* heap, FILE *stream) {
  fputc('[', stream);
  for (size_t i = 0; i < heap->size; ++i) {
    lp_polynomial_print(heap->data[i], stream);
    if (i < heap->size - 1)
      fputc(',', stream);
  }
  fputc(']', stream);
}