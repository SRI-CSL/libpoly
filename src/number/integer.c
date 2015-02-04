/*
 * coefficient.c
 *
 *  Created on: Nov 5, 2013
 *      Author: dejan
 */

#include "number/integer.h"

int integer_print_matrix(const lp_integer_t* c, size_t m, size_t n, FILE* out) {
  size_t i, j;
  int len = 0;
  for (i = 0; i < m; ++ i) {
    for (j = 0; j < n; ++ j) {
      len += gmp_fprintf(out, "%4Zd", c + i*m + j);
    }
    len += fprintf(out, "\n");
  }
  return len;
}

static
lp_int_ring integer_ring_create(const lp_integer_t* M, int is_prime) {

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

static
void integer_ring_destroy(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  mpz_clear(&nonconst->M);
  mpz_clear(&nonconst->ub);
  mpz_clear(&nonconst->lb);
  free(nonconst);
}

static
void integer_ring_attach(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    ++ nonconst->ref_count;
  }
}

static
void integer_ring_detach(lp_int_ring K) {
  lp_int_ring_t* nonconst = (lp_int_ring_t*) K;
  if (nonconst) {
    assert(nonconst->ref_count > 0);
    if (--nonconst->ref_count == 0) {
      integer_ring_destroy(nonconst);
    }
  }
}

static
int integer_ring_equal(lp_int_ring K1, lp_int_ring K2) {
  if (K1 == K2) {
    return 1;
  } else if (K1 && K2) {
    return mpz_cmp(&K1->M, &K2->M) == 0;
  } else {
    return 0;
  }
}


static
int integer_ring_print(lp_int_ring K, FILE* out) {
  int len = 0;
  len += fprintf(out, "Z");
  if (K) {
    len += fprintf(out, " mod ");
    len += integer_print(&K->M, out);
  }
  return len;
}

static
char* integer_ring_to_string(lp_int_ring K) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  integer_ring_print(K, f);
  fclose(f);
  return str;
}

lp_int_ring_ops_t lp_int_ring_ops = {
    integer_ring_create,
    integer_ring_attach,
    integer_ring_detach,
    integer_ring_equal,
    integer_ring_print,
    integer_ring_to_string
};

lp_int_ring lp_Z = 0;

static
void integer_construct_from_string(lp_int_ring K, lp_integer_t* c, const char* x, int base) {
  mpz_init_set_str(c, x, base);
  integer_ring_normalize(K, c);
}


static
char* integer_to_string(const lp_integer_t* c) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  integer_print(c, f);
  fclose(f);
  return str;
}

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

