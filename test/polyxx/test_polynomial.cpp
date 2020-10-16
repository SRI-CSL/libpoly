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
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  Polynomial r = resultant(p, q);
  CHECK(r == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
}

TEST_CASE("polynomial::discriminant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial d = discriminant(p);
  CHECK(d == 11337408 * pow(y, 5) - 35095680 * pow(y, 4) + 34197120 * pow(y, 3) - 14999040 * pow(y, 2) + 3099840 * y - 246656);
}
TEST_CASE("polynomial::psc") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  auto tmp = psc(p, q);

  CHECK(tmp.size() == 6);
  CHECK(tmp[0] == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
  CHECK(tmp[1] == 416745 * y * y * y - 416745 * y * y + 138915 * y - 15435);
  CHECK(tmp[2] == Integer(0));
  CHECK(tmp[3] == Integer(0));
  CHECK(tmp[4] == Integer(-45));
  CHECK(tmp[5] == Integer(7));
}