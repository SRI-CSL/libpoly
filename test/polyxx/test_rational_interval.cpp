#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("rational_interval::constructors") {
  RationalInterval ri1;
  CHECK(is_point(ri1));
  CHECK(get_point(ri1) == Rational(0));
  CHECK(get_lower(ri1) == Rational(0));
  CHECK(get_upper(ri1) == Rational(0));

  RationalInterval ri2(Rational(1));
  CHECK(is_point(ri2));
  CHECK(get_point(ri2) == Rational(1));
  CHECK(get_lower(ri2) == Rational(1));
  CHECK(get_upper(ri2) == Rational(1));

  RationalInterval ri3(Rational(1), Rational(2));
  CHECK_FALSE(is_point(ri3));
  CHECK(get_lower(ri3) == Rational(1));
  CHECK(get_upper(ri3) == Rational(2));
}

TEST_CASE("rational_interval::contains_zero") {
  CHECK_FALSE(contains_zero(RationalInterval(-2, -1)));
  CHECK(contains_zero(RationalInterval(-1, 1)));
  CHECK_FALSE(contains_zero(RationalInterval(1, 2)));
}

TEST_CASE("rational_interval::contains") {
  CHECK_FALSE(contains(RationalInterval(1, 3), Rational(0)));
  CHECK(contains(RationalInterval(1, 3), Rational(2)));
  CHECK_FALSE(contains(RationalInterval(1, 3), Rational(4)));

  CHECK_FALSE(contains(RationalInterval(1, 3), Value(Rational(0))));
  CHECK(contains(RationalInterval(1, 3), Value(Rational(2))));
  CHECK_FALSE(contains(RationalInterval(1, 3), Value(Rational(4))));
}

TEST_CASE("rational_interval::sgn") {
  CHECK(sgn(RationalInterval(-3, -2)) == -1);
  CHECK(sgn(RationalInterval(-3, 0)) == -1);
  CHECK(sgn(RationalInterval(-3, 2)) == 0);
  CHECK(sgn(RationalInterval(0, 2)) == 1);
  CHECK(sgn(RationalInterval(1, 2)) == 1);
}

TEST_CASE("rational_interval::operator<<") {
  RationalInterval ri(-3, 2);
  std::stringstream out;
  out << ri;
  CHECK(out.str() == "( -3 ; 2 )");
}
