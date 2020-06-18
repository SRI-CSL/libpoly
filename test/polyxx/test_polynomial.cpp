#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("polynomial::properties") {
  Variable x("x");
  Variable y("y");
  Polynomial p1;
  Polynomial p2(Integer(5));
  Polynomial p3 = 3*x;
  Polynomial p4 = 3*x*y*y + 7*x;

  CHECK(is_zero(p1));
  CHECK_FALSE(is_zero(p2));
  CHECK_FALSE(is_zero(p3));
  CHECK_FALSE(is_zero(p4));

  CHECK(is_constant(p1));
  CHECK(is_constant(p2));
  CHECK_FALSE(is_constant(p3));
  CHECK_FALSE(is_constant(p4));

  CHECK_FALSE(is_linear(p1));
  CHECK_FALSE(is_linear(p2));
  CHECK(is_linear(p3));
  CHECK_FALSE(is_linear(p4));

  CHECK(is_lc_constant(p1));
  CHECK(is_lc_constant(p2));
  CHECK(is_lc_constant(p3));
  CHECK_FALSE(is_lc_constant(p4));

  CHECK(degree(p1) == 0);
  CHECK(degree(p2) == 0);
  CHECK(degree(p3) == 1);
  CHECK(degree(p4) == 2);

  CHECK(is_univariate(p3));
}

TEST_CASE("polynomial::resultant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 7*pow(x, 6) + 14*pow(x, 5) + 7*y - 1;
  Polynomial q = 42*pow(x, 5) + 70*pow(x, 4);
  Polynomial r = resultant(p, q);
  CHECK(r == 92254156521408*pow(y,5) - 461361174686720*pow(y,4) + 244807578081920*pow(y,3) - 51113953954560*pow(y,2) + 4803956911040*y - 170197631744);
}

TEST_CASE("polynomial::discriminant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 7*pow(x, 6) + 14*pow(x, 5) + 7*y - 1;
  Polynomial d = discriminant(p);
  CHECK(d == 13179165217344*pow(y,5) - 65908739240960*pow(y,4) + 34972511154560*pow(y,3) - 7301993422080*pow(y,2) + 686279558720*y - 24313947392);
}
TEST_CASE("polynomial::psc") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 7*pow(x, 6) + 14*pow(x, 5) + 7*y - 1;
  Polynomial q = 42*pow(x, 5) + 70*pow(x, 4);
  auto tmp = psc(p, q);

  CHECK(tmp.size() == 6);
  CHECK(tmp[0] == 92254156521408*pow(y,5) - 461361174686720*pow(y,4) + 244807578081920*pow(y,3) - 51113953954560*pow(y,2) + 4803956911040*y - 170197631744);
  CHECK(tmp[1] == 174327582240*y*y*y-74711820960*y*y+10673117280*y-508243680);
  CHECK(tmp[2] == Integer(0));
  CHECK(tmp[3] == Integer(0));
  CHECK(tmp[4] == Integer(-6860));
  CHECK(tmp[5] == Integer(42));
}