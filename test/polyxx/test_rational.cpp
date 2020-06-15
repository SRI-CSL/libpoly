#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("rational::constructors") {
  CHECK(Rational() == Rational());
  CHECK(Rational() == Rational(0, 1));
  CHECK(Rational(13) == Rational(26, 2));
  CHECK(Rational(Integer(7)) == Rational(21, 3));
  CHECK(Rational(3.0) == Rational(6, 2));
}

TEST_CASE("rational::to_double") {
  CHECK(to_double(Rational(1, 2)) == 0.5);
  CHECK(to_double(Rational(1)) == 1.0);
  CHECK(to_double(Rational(2)) == 2.0);
}

TEST_CASE("rational::sgn") {
  CHECK(sgn(Rational(-10)) == -1);
  CHECK(sgn(Rational(-1)) == -1);
  CHECK(sgn(Rational()) == 0);
  CHECK(sgn(Rational(1)) == 1);
  CHECK(sgn(Rational(10)) == 1);
}

TEST_CASE("rational::operator==") {
  CHECK_FALSE(Rational(1) == Rational(2));
  CHECK_FALSE(Rational(1) == Integer(2));
  CHECK_FALSE(Integer(1) == Rational(2));
  CHECK(Rational(1) == Rational(1));
  CHECK(Rational(1) == Integer(1));
  CHECK(Integer(1) == Rational(1));
  CHECK_FALSE(Rational(2) == Rational(1));
  CHECK_FALSE(Rational(2) == Integer(1));
  CHECK_FALSE(Integer(2) == Rational(1));
}

TEST_CASE("rational::operator!=") {
  CHECK(Rational(1) != Rational(2));
  CHECK(Rational(1) != Integer(2));
  CHECK(Integer(1) != Rational(2));
  CHECK_FALSE(Rational(1) != Rational(1));
  CHECK_FALSE(Rational(1) != Integer(1));
  CHECK_FALSE(Integer(1) != Rational(1));
  CHECK(Rational(2) != Rational(1));
  CHECK(Rational(2) != Integer(1));
  CHECK(Integer(2) != Rational(1));
}

TEST_CASE("rational::operator<") {
  CHECK(Rational(1) < Rational(2));
  CHECK(Rational(1) < Integer(2));
  CHECK(Integer(1) < Rational(2));
  CHECK_FALSE(Rational(1) < Rational(1));
  CHECK_FALSE(Rational(1) < Integer(1));
  CHECK_FALSE(Integer(1) < Rational(1));
  CHECK_FALSE(Rational(2) < Rational(1));
  CHECK_FALSE(Rational(2) < Integer(1));
  CHECK_FALSE(Integer(2) < Rational(1));
}

TEST_CASE("rational::operator<=") {
  CHECK(Rational(1) <= Rational(2));
  CHECK(Rational(1) <= Integer(2));
  CHECK(Integer(1) <= Rational(2));
  CHECK(Rational(1) <= Rational(1));
  CHECK(Rational(1) <= Integer(1));
  CHECK(Integer(1) <= Rational(1));
  CHECK_FALSE(Rational(2) <= Rational(1));
  CHECK_FALSE(Rational(2) <= Integer(1));
  CHECK_FALSE(Integer(2) <= Rational(1));
}

TEST_CASE("rational::operator>") {
  CHECK_FALSE(Rational(1) > Rational(2));
  CHECK_FALSE(Rational(1) > Integer(2));
  CHECK_FALSE(Integer(1) > Rational(2));
  CHECK_FALSE(Rational(1) > Rational(1));
  CHECK_FALSE(Rational(1) > Integer(1));
  CHECK_FALSE(Integer(1) > Rational(1));
  CHECK(Rational(2) > Rational(1));
  CHECK(Rational(2) > Integer(1));
  CHECK(Integer(2) > Rational(1));
}

TEST_CASE("rational::operator>=") {
  CHECK_FALSE(Rational(1) >= Rational(2));
  CHECK_FALSE(Rational(1) >= Integer(2));
  CHECK_FALSE(Integer(1) >= Rational(2));
  CHECK(Rational(1) >= Rational(1));
  CHECK(Rational(1) >= Integer(1));
  CHECK(Integer(1) >= Rational(1));
  CHECK(Rational(2) >= Rational(1));
  CHECK(Rational(2) >= Integer(1));
  CHECK(Integer(2) >= Rational(1));
}

TEST_CASE("rational::swap") {
  Rational a(1);
  Rational b(2);
  swap(a, b);
  CHECK(a == Rational(2));
  CHECK(b == Rational(1));
}

TEST_CASE("rational::operator+") {
  CHECK(Rational(1) + Rational(2) == Rational(3));
  CHECK(Rational(1) + Integer(2) == Rational(3));
  CHECK(Integer(1) + Rational(2) == Rational(3));
}
TEST_CASE("rational::operator-") {
  CHECK(Rational(1) - Rational(2) == Rational(-1));
  CHECK(-Rational(1) == Rational(-1));
}

TEST_CASE("rational::inverse") {
  CHECK(inverse(Rational(1)) == Rational(1));
  CHECK(inverse(Rational(Integer(1), Integer(2))) == Rational(2));
  CHECK(inverse(Rational(2)) == Rational(Integer(1), Integer(2)));
}
TEST_CASE("rational::operator*") {
  CHECK(Rational(1) * Rational(2) == Rational(2));
}
TEST_CASE("rational::mul_2exp") {
  CHECK(mul_2exp(Rational(2), 2) == Rational(8));
}
TEST_CASE("rational::pow") { CHECK(pow(Rational(2), 3) == Rational(8)); }
TEST_CASE("rational::operator/") {
  CHECK(Rational(1) / Rational(1) == Rational(1));
  CHECK(Rational(1) / Rational(2) == Rational(Integer(1), Integer(2)));
}

TEST_CASE("rational::div_2exp") {
  CHECK(div_2exp(Rational(2), 2) == Rational(Integer(1), Integer(2)));
}

TEST_CASE("rational::numerator") {
  CHECK(numerator(Rational(1)) == Rational(1));
  CHECK(numerator(Rational(Integer(1), Integer(2))) == Rational(1));
  CHECK(numerator(Rational(2)) == Rational(2));
}
TEST_CASE("rational::denominator") {
  CHECK(denominator(Rational(1)) == Rational(1));
  CHECK(denominator(Rational(Integer(1), Integer(2))) == Rational(2));
  CHECK(denominator(Rational(2)) == Rational(1));
}

TEST_CASE("rational::is_integer") {
  CHECK(is_integer(Rational(1)));
  CHECK_FALSE(is_integer(Rational(Integer(1), Integer(2))));
  CHECK(is_integer(Rational(2)));
}

TEST_CASE("rational::ceil") {
  CHECK(ceil(Rational(1)) == Rational(1));
  CHECK(ceil(Rational(Integer(1), Integer(2))) == Rational(1));
}
TEST_CASE("rational::floor") {
  CHECK(floor(Rational(1)) == Rational(1));
  CHECK(floor(Rational(Integer(1), Integer(2))) == Rational());
}
