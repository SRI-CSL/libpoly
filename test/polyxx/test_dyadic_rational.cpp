#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("dyadic_rational::constructors") {
  CHECK(DyadicRational() == DyadicRational());
  CHECK(DyadicRational() == DyadicRational(0, 1));
  CHECK(DyadicRational(0.0) == DyadicRational(0, 1));
  CHECK(DyadicRational(0.5) == DyadicRational(1, 1));
  CHECK(DyadicRational(2) == DyadicRational(Integer(2)));
}

TEST_CASE("dyadic_rational::rational") {
  CHECK(Rational(DyadicRational(2)) == Rational(Integer(2)));
}

TEST_CASE("dyadic_rational::to_double") {
  CHECK(to_double(DyadicRational(Integer(1))) == 1.0);
  CHECK(to_double(DyadicRational(Integer(2))) == 2.0);
  CHECK(to_double(DyadicRational(1, 2)) == 0.25);
  CHECK(to_double(DyadicRational(3, 2)) == 0.75);
}

TEST_CASE("dyadic_rational::sgn") {
  CHECK(sgn(DyadicRational(-10)) == -1);
  CHECK(sgn(DyadicRational(-1)) == -1);
  CHECK(sgn(DyadicRational(0, 1)) == 0);
  CHECK(sgn(DyadicRational(1)) == 1);
  CHECK(sgn(DyadicRational(10)) == 1);
}

TEST_CASE("dyadic_rational::operator==") {
  CHECK_FALSE(DyadicRational(1) == DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) == Integer(2));
  CHECK_FALSE(Integer(1) == DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) == Rational(2));
  CHECK_FALSE(Rational(1) == DyadicRational(2));
  CHECK(DyadicRational(1) == DyadicRational(1));
  CHECK(DyadicRational(1) == Integer(1));
  CHECK(Integer(1) == DyadicRational(1));
  CHECK(DyadicRational(1) == Rational(1));
  CHECK(Rational(1) == DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) == DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) == Integer(1));
  CHECK_FALSE(Integer(2) == DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) == Rational(1));
  CHECK_FALSE(Rational(2) == DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator!=") {
  CHECK(DyadicRational(1) != DyadicRational(2));
  CHECK(DyadicRational(1) != Integer(2));
  CHECK(Integer(1) != DyadicRational(2));
  CHECK(DyadicRational(1) != Rational(2));
  CHECK(Rational(1) != DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) != DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) != Integer(1));
  CHECK_FALSE(Integer(1) != DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) != Rational(1));
  CHECK_FALSE(Rational(1) != DyadicRational(1));
  CHECK(DyadicRational(2) != DyadicRational(1));
  CHECK(DyadicRational(2) != Integer(1));
  CHECK(Integer(2) != DyadicRational(1));
  CHECK(DyadicRational(2) != Rational(1));
  CHECK(Rational(2) != DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator<") {
  CHECK(DyadicRational(1) < DyadicRational(2));
  CHECK(DyadicRational(1) < Integer(2));
  CHECK(Integer(1) < DyadicRational(2));
  CHECK(DyadicRational(1) < Rational(2));
  CHECK(Rational(1) < DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) < DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) < Integer(1));
  CHECK_FALSE(Integer(1) < DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) < Rational(1));
  CHECK_FALSE(Rational(1) < DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) < DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) < Integer(1));
  CHECK_FALSE(Integer(2) < DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) < Rational(1));
  CHECK_FALSE(Rational(2) < DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator<=") {
  CHECK(DyadicRational(1) <= DyadicRational(2));
  CHECK(DyadicRational(1) <= Integer(2));
  CHECK(Integer(1) <= DyadicRational(2));
  CHECK(DyadicRational(1) <= Rational(2));
  CHECK(Rational(1) <= DyadicRational(2));
  CHECK(DyadicRational(1) <= DyadicRational(1));
  CHECK(DyadicRational(1) <= Integer(1));
  CHECK(Integer(1) <= DyadicRational(1));
  CHECK(DyadicRational(1) <= Rational(1));
  CHECK(Rational(1) <= DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) <= DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) <= Integer(1));
  CHECK_FALSE(Integer(2) <= DyadicRational(1));
  CHECK_FALSE(DyadicRational(2) <= Rational(1));
  CHECK_FALSE(Rational(2) <= DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator>") {
  CHECK_FALSE(DyadicRational(1) > DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) > Integer(2));
  CHECK_FALSE(Integer(1) > DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) > Rational(2));
  CHECK_FALSE(Rational(1) > DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) > DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) > Integer(1));
  CHECK_FALSE(Integer(1) > DyadicRational(1));
  CHECK_FALSE(DyadicRational(1) > Rational(1));
  CHECK_FALSE(Rational(1) > DyadicRational(1));
  CHECK(DyadicRational(2) > DyadicRational(1));
  CHECK(DyadicRational(2) > Integer(1));
  CHECK(Integer(2) > DyadicRational(1));
  CHECK(DyadicRational(2) > Rational(1));
  CHECK(Rational(2) > DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator>=") {
  CHECK_FALSE(DyadicRational(1) >= DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) >= Integer(2));
  CHECK_FALSE(Integer(1) >= DyadicRational(2));
  CHECK_FALSE(DyadicRational(1) >= Rational(2));
  CHECK_FALSE(Rational(1) >= DyadicRational(2));
  CHECK(DyadicRational(1) >= DyadicRational(1));
  CHECK(DyadicRational(1) >= Integer(1));
  CHECK(Integer(1) >= DyadicRational(1));
  CHECK(DyadicRational(1) >= Rational(1));
  CHECK(Rational(1) >= DyadicRational(1));
  CHECK(DyadicRational(2) >= DyadicRational(1));
  CHECK(DyadicRational(2) >= Integer(1));
  CHECK(Integer(2) >= DyadicRational(1));
  CHECK(DyadicRational(2) >= Rational(1));
  CHECK(Rational(2) >= DyadicRational(1));
}

TEST_CASE("dyadic_rational::swap") {
  DyadicRational a(1);
  DyadicRational b(2);
  swap(a, b);
  CHECK(a == DyadicRational(2));
  CHECK(b == DyadicRational(1));
}

TEST_CASE("dyadic_rational::operator+") {
  CHECK(DyadicRational(1) + DyadicRational(2) == DyadicRational(3));
  CHECK(DyadicRational(1) + Integer(2) == DyadicRational(3));
  CHECK(Integer(1) + DyadicRational(2) == DyadicRational(3));
}
TEST_CASE("dyadic_rational::operator-") {
  CHECK(DyadicRational(1) - DyadicRational(2) == DyadicRational(-1));
  CHECK(-DyadicRational(1) == DyadicRational(-1));
}

TEST_CASE("dyadic_rational::operator*") {
  CHECK(DyadicRational(1) * DyadicRational(2) == DyadicRational(2));
}
TEST_CASE("dyadic_rational::mul_2exp") {
  CHECK(mul_2exp(DyadicRational(2), 2) == DyadicRational(8));
}
TEST_CASE("dyadic_rational::pow") {
  CHECK(pow(DyadicRational(2), 3) == DyadicRational(8));
}

TEST_CASE("dyadic_rational::div_2exp") {
  CHECK(div_2exp(DyadicRational(2), 1) == DyadicRational(Integer(1)));
}

TEST_CASE("dyadic_rational::numerator") {
  CHECK(numerator(DyadicRational(1)) == Integer(1));
  CHECK(numerator(DyadicRational(3, 4)) == Integer(3));
  CHECK(numerator(DyadicRational(2, 4)) == Integer(1));
}
TEST_CASE("dyadic_rational::denominator") {
  CHECK(denominator(DyadicRational(1)) == Integer(1));
  CHECK(denominator(DyadicRational(1, 2)) == Integer(4));
}

TEST_CASE("dyadic_rational::is_integer") {
  CHECK(is_integer(DyadicRational(1)));
  CHECK(is_integer(DyadicRational(2)));
  CHECK_FALSE(is_integer(DyadicRational(2, 4)));
  CHECK(is_integer(DyadicRational(4, 2)));
  CHECK_FALSE(is_integer(DyadicRational(3, 4)));
}

TEST_CASE("dyadic_rational::ceil") {
  CHECK(ceil(DyadicRational(1)) == Integer(1));
  CHECK(ceil(DyadicRational(1, 3)) == Integer(1));
  CHECK(ceil(DyadicRational(7, 2)) == Integer(2));
}
TEST_CASE("dyadic_rational::floor") {
  CHECK(floor(DyadicRational(1)) == Integer(1));
  CHECK(floor(DyadicRational(1, 3)) == Integer(0));
  CHECK(floor(DyadicRational(7, 2)) == Integer(1));
}

TEST_CASE("dyadic_rational::operator<<") {
  DyadicRational dr(1, 3);
  std::stringstream out;
  out << dr;
  CHECK(out.str() == "1/8");
}
