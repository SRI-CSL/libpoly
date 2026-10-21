#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../polyxx/doctest.h"

#include <gmp.h>
#include <poly.h>
#include <integer.h>
#include <upolynomial.h>

TEST_CASE("upolynomial::gcd_shift_past_unsigned") {
  // g = x^64 - 2^(2^26): the heuristic would shift by degree * n > 2^32
  // bits, past an unsigned, so the gcd comes from the subresultant algorithm
  lp_integer_t g_coeffs[65];
  for (int i = 0; i < 65; i ++) {
    lp_integer_construct_from_int(lp_Z, &g_coeffs[i], 0);
  }
  lp_integer_assign_int(lp_Z, &g_coeffs[0], 1);
  mpz_mul_2exp(&g_coeffs[0], &g_coeffs[0], 1 << 26);
  mpz_neg(&g_coeffs[0], &g_coeffs[0]);
  lp_integer_assign_int(lp_Z, &g_coeffs[64], 1);
  lp_upolynomial_t* g = lp_upolynomial_construct(lp_Z, 64, g_coeffs);

  int p_coeffs[2] = { 1, 1 };
  int q_coeffs[2] = { 2, 1 };
  lp_upolynomial_t* p = lp_upolynomial_construct_from_int(lp_Z, 1, p_coeffs);
  lp_upolynomial_t* q = lp_upolynomial_construct_from_int(lp_Z, 1, q_coeffs);
  lp_upolynomial_t* gp = lp_upolynomial_mul(g, p);
  lp_upolynomial_t* gq = lp_upolynomial_mul(g, q);

  lp_upolynomial_t* gcd = lp_upolynomial_gcd(gp, gq);
  CHECK(lp_upolynomial_cmp(gcd, g) == 0);

  lp_upolynomial_delete(gcd);
  lp_upolynomial_delete(gq);
  lp_upolynomial_delete(gp);
  lp_upolynomial_delete(q);
  lp_upolynomial_delete(p);
  lp_upolynomial_delete(g);
  for (int i = 0; i < 65; i ++) {
    lp_integer_destruct(&g_coeffs[i]);
  }
}
