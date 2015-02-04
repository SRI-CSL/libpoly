/*
 * coefficient.c
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#include "number/integer.h"

lp_int_ring lp_integer_ring_create(const lp_integer_t* M, int is_prime) {

  assert(mpz_sgn(M) > 0);

  lp_integer_t tmp;
  mpz_init(&tmp);

  lp_int_ring_t* K = malloc(sizeof(lp_int_ring_t));

  K->is_prime = is_prime;
  assert(!!is_prime == !!mpz_probab_prime_p(M, 25));
  K->ref_count = 1;

  mpz_init_set(&K->M, M);
  mpz_init(&K->ub);
  mpz_tdiv_q_2exp(&K->ub, M, 1);

  mpz_init(&K->lb);
  mpz_sub_ui(&K->lb, M, 1); // M-1
  mpz_tdiv_q_2exp(&tmp, &K->lb, 1); // [M-1/2]
  mpz_neg(&K->lb, &tmp);

  mpz_clear(&tmp);

  return K;
}

void lp_integer_ring_destroy(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  mpz_clear(&nonconst->M);
  mpz_clear(&nonconst->ub);
  mpz_clear(&nonconst->lb);
  free(nonconst);
}

void lp_integer_ring_attach(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    ++ nonconst->ref_count;
  }
}

void lp_integer_ring_detach(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    assert(nonconst->ref_count > 0);
    if (--nonconst->ref_count == 0) {
      lp_integer_ring_destroy(nonconst);
    }
  }
}

int lp_integer_ring_equal(lp_int_ring K1, lp_int_ring K2) {
  if (K1 == K2) {
    return 1;
  } else if (K1 && K2) {
    return mpz_cmp(&K1->M, &K2->M) == 0;
  } else {
    return 0;
  }
}


int lp_integer_ring_print(lp_int_ring K, FILE* out) {
  int len = 0;
  len += fprintf(out, "Z");
  if (K) {
    len += fprintf(out, " mod ");
    len += integer_print(&K->M, out);
  }
  return len;
}

char* lp_integer_ring_to_string(lp_int_ring K) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_integer_ring_print(K, f);
  fclose(f);
  return str;
}

lp_int_ring lp_Z = 0;

const lp_integer_ops_struct lp_integer_ops = {
    integer_construct_from_int,
    integer_construct_from_string,
    integer_construct_copy,
    integer_assign,
    integer_assign_int,
    integer_destruct,
    integer_print,
    integer_bits,
    integer_print_matrix,
    integer_to_string,
    integer_to_int,
    integer_is_prime,
    integer_is_zero,
    integer_in_ring,
    integer_sgn,
    integer_cmp,
    integer_cmp_int,
    integer_divides,
    integer_swap,
    integer_inc,
    integer_dec,
    integer_add,
    integer_sub,
    integer_neg,
    integer_abs,
    integer_inv,
    integer_mul,
    integer_mul_int,
    integer_mul_pow2,
    integer_pow,
    integer_sqrt_Z,
    integer_add_mul,
    integer_add_mul_int,
    integer_sub_mul,
    integer_div_exact,
    integer_div_Z,
    integer_rem_Z,
    integer_div_rem_Z,
    integer_div_rem_pow2_Z,
    integer_gcd_Z,
    integer_lcm_Z
};

