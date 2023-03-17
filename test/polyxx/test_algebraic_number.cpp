#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("algebraic_number::move") {
  AlgebraicNumber a(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1));
  AlgebraicNumber a2(std::move(a));
}

TEST_CASE("algebraic_number::swap") {
  AlgebraicNumber a(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1));
  AlgebraicNumber b(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2));
  swap(a, b);
  CHECK(a == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2)));
  CHECK(b == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1)));
}

TEST_CASE("algebraic_number::sgn") {
  CHECK(sgn(AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1))) ==
        -1);
  CHECK(sgn(AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2))) ==
        1);
}

TEST_CASE("algebraic_number::to_double") {
  double d1 = to_double(
      AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1)));
  CHECK(-1.5 < d1);
  CHECK(d1 < -1.4);
  double d2 =
      to_double(AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2)));
  CHECK(1.4 < d2);
  CHECK(d2 < 1.5);
}
TEST_CASE("algebraic_number::to_rational_approximation") {
  Rational r1 = to_rational_approximation(
      AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1)));
  CHECK(Rational(-2) < r1);
  CHECK(r1 < Rational(-1));
  Rational r2 = to_rational_approximation(
      AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2)));
  CHECK(Rational(1) < r2);
  CHECK(r2 < Rational(2));
}

TEST_CASE("algebraic_number::is_rational") {
  CHECK_FALSE(is_rational(
      AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1))));
  CHECK(
      is_rational(AlgebraicNumber(UPolynomial({1, 3}), DyadicInterval(-2, 0))));
  CHECK(is_rational(
      AlgebraicNumber(UPolynomial({3, 1}), DyadicInterval(-4, -2))));
}
TEST_CASE("algebraic_number::is_integer") {
  CHECK_FALSE(is_integer(
      AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1))));
  CHECK_FALSE(
      is_integer(AlgebraicNumber(UPolynomial({1, 3}), DyadicInterval(-2, 0))));
  CHECK(
      is_integer(AlgebraicNumber(UPolynomial({3, 1}), DyadicInterval(-4, -2))));
}

TEST_CASE("algebraic_number::to_rational") {
  {
    AlgebraicNumber an(UPolynomial({1, 3}), DyadicInterval(-2, 0));
    CHECK(is_rational(an));
    CHECK(to_rational_approximation(an) == Rational(-1, 3));
  }
  {
    AlgebraicNumber an(UPolynomial({3, 1}), DyadicInterval(-4, -2));
    CHECK(is_rational(an));
    CHECK(to_rational_approximation(an) == Rational(-3));
  }
}

TEST_CASE("algebraic_number::ceil") {
  CHECK(ceil(AlgebraicNumber(UPolynomial({-2, 0, 1}),
                             DyadicInterval(-2, -1))) == Integer(-1));
  CHECK(ceil(AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2))) ==
        Integer(2));
}
TEST_CASE("algebraic_number::floor") {
  CHECK(floor(AlgebraicNumber(UPolynomial({-2, 0, 1}),
                              DyadicInterval(-2, -1))) == Integer(-2));
  CHECK(floor(AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2))) ==
        Integer(1));
}
TEST_CASE("algebraic_number::neg") {
  AlgebraicNumber two(UPolynomial({-2, 1}), DyadicInterval(1, 3));
  CHECK(-two == Rational(-2));

  AlgebraicNumber sqrt2(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2));
  CHECK(-(sqrt2*sqrt2) == Rational(-2));

  AlgebraicNumber msqrt2(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1));
  CHECK(-sqrt2 == msqrt2);
  CHECK(sqrt2 == -msqrt2);
}
TEST_CASE("algebraic_number::operator/") {
  AlgebraicNumber sqrt2(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2));
  CHECK(is_one(sqrt2 / sqrt2));
  CHECK(sqrt2 / -sqrt2 == Rational(-1));
  CHECK(-(sqrt2 / -sqrt2) == Rational(1));
}

TEST_CASE("algebraic_number::operator<<") {
  AlgebraicNumber a(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1));
  std::stringstream out;
  out << a;
  CHECK(out.str() == "<1*x^2 + (-2), (-3/2, -5/4)>");
}
