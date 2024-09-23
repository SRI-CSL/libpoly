#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("interval::constructors") {
  Interval zero;
  CHECK(is_point(zero));
  CHECK(get_point(zero) == Value(Integer(0)));
  CHECK(get_lower(zero) == Value(Integer(0)));
  CHECK(get_upper(zero) == Value(Integer(0)));

  Interval point(Value(Integer(3)));
  CHECK(is_point(point));
  CHECK(get_point(point) == Value(Integer(3)));
  CHECK(get_lower(point) == Value(Integer(3)));
  CHECK(get_upper(point) == Value(Integer(3)));

  Interval open(3, 7);
  CHECK_FALSE(is_point(open));
  CHECK(get_lower(open) == Value(Integer(3)));
  CHECK(get_upper(open) == Value(Integer(7)));

  Interval closed(3, 7);
  CHECK_FALSE(is_point(closed));
  CHECK(get_lower(closed) == Value(Integer(3)));
  CHECK(get_upper(closed) == Value(Integer(7)));

  closed.set_lower(4, true);
  CHECK_FALSE(is_point(closed));
  CHECK(get_lower(closed) == Value(Integer(4)));
  CHECK(get_upper(closed) == Value(Integer(7)));

  closed.set_upper(5, false);
  CHECK_FALSE(is_point(closed));
  CHECK(get_lower(closed) == Value(Integer(4)));
  CHECK(get_upper(closed) == Value(Integer(5)));

  closed.set_lower(5, false);
  CHECK(is_point(closed));
  CHECK(get_point(closed) == Value(Integer(5)));
  CHECK(get_lower(closed) == Value(Integer(5)));
  CHECK(get_upper(closed) == Value(Integer(5)));
}

TEST_CASE("interval::contains") {
  Interval i(1,3);
  CHECK(contains(i, Value(AlgebraicNumber(2))));
  CHECK(contains(i, Value(DyadicRational(5,2))));
  CHECK(contains(i, Value(Integer(2))));
  CHECK(contains(i, Value(Rational(5,2))));
}

TEST_CASE("interval::contains_int") {
  CHECK_FALSE(contains_int(Interval(1,2)));
  CHECK(contains_int(Interval(1,3)));
  CHECK(contains_int(Interval()));
  CHECK(contains_int(Interval(1)));
  CHECK_FALSE(contains_int(Interval(Value(Rational(3,2)))));
  CHECK(contains_int(Interval(1, false, 2, true)));
  CHECK(contains_int(Interval(1, true, 2, false)));
  CHECK_FALSE(contains_int(Interval(Rational(1,4), Rational(3,4))));
}

TEST_CASE("interval::count_int") {
  CHECK(count_int(Interval(0,1)) == 0);
  CHECK(count_int(Interval(0,2)) == 1);
  CHECK(count_int(Interval()) == 1);
  CHECK(count_int(Interval(Rational(1,2))) == 0);
  CHECK(count_int(Interval(Rational(3,4), Rational(5,4))) == 1);
  CHECK(count_int(Interval(0, false, 2, false)) == 3);
  CHECK(count_int(Interval(0,LONG_MAX)) == LONG_MAX - 1);
  CHECK(count_int(Interval(0,false, LONG_MAX, true)) == LONG_MAX);
  CHECK(count_int(Interval(0,false, LONG_MAX, false)) == LONG_MAX);
  CHECK(count_int(Interval(-1,LONG_MAX)) == LONG_MAX);
  CHECK(count_int(Interval(-2,LONG_MAX)) == LONG_MAX);
  CHECK(count_int(Interval(LONG_MIN, LONG_MAX)) == LONG_MAX);
}

TEST_CASE("interval::pick_value") {
  Interval i(1,3);
  CHECK(contains(i, pick_value(i)));
}

TEST_CASE("interval::operator==") {
  CHECK_FALSE(Interval(1,2) == Interval(3,6));
  CHECK_FALSE(Interval(1,3) == Interval(3,6));
  CHECK_FALSE(Interval(1,4) == Interval(3,6));
  CHECK_FALSE(Interval(1,6) == Interval(3,6));
  CHECK_FALSE(Interval(1,8) == Interval(3,6));
  CHECK_FALSE(Interval(3,5) == Interval(3,6));
  CHECK(Interval(3,6) == Interval(3,6));
  CHECK_FALSE(Interval(3,8) == Interval(3,6));
  CHECK_FALSE(Interval(4,5) == Interval(3,6));
  CHECK_FALSE(Interval(4,6) == Interval(3,6));
  CHECK_FALSE(Interval(4,8) == Interval(3,6));
  CHECK_FALSE(Interval(6,8) == Interval(3,6));
  CHECK_FALSE(Interval(7,8) == Interval(3,6));
}

TEST_CASE("interval::operator!=") {
  CHECK(Interval(1,2) != Interval(3,6));
  CHECK(Interval(1,3) != Interval(3,6));
  CHECK(Interval(1,4) != Interval(3,6));
  CHECK(Interval(1,6) != Interval(3,6));
  CHECK(Interval(1,8) != Interval(3,6));
  CHECK(Interval(3,5) != Interval(3,6));
  CHECK_FALSE(Interval(3,6) != Interval(3,6));
  CHECK(Interval(3,8) != Interval(3,6));
  CHECK(Interval(4,5) != Interval(3,6));
  CHECK(Interval(4,6) != Interval(3,6));
  CHECK(Interval(4,8) != Interval(3,6));
  CHECK(Interval(6,8) != Interval(3,6));
  CHECK(Interval(7,8) != Interval(3,6));
}

TEST_CASE("interval::operator<") {
  CHECK(Interval(1,2) < Interval(3,6));
  CHECK(Interval(1,3) < Interval(3,6));
  CHECK(Interval(1,4) < Interval(3,6));
  CHECK(Interval(1,6) < Interval(3,6));
  CHECK(Interval(1,8) < Interval(3,6));
  CHECK(Interval(3,5) < Interval(3,6));
  CHECK_FALSE(Interval(3,6) < Interval(3,6));
  CHECK_FALSE(Interval(3,8) < Interval(3,6));
  CHECK_FALSE(Interval(4,5) < Interval(3,6));
  CHECK_FALSE(Interval(4,6) < Interval(3,6));
  CHECK_FALSE(Interval(4,8) < Interval(3,6));
  CHECK_FALSE(Interval(6,8) < Interval(3,6));
  CHECK_FALSE(Interval(7,8) < Interval(3,6));
}

TEST_CASE("interval::operator<=") {
  CHECK(Interval(1,2) <= Interval(3,6));
  CHECK(Interval(1,3) <= Interval(3,6));
  CHECK(Interval(1,4) <= Interval(3,6));
  CHECK(Interval(1,6) <= Interval(3,6));
  CHECK(Interval(1,8) <= Interval(3,6));
  CHECK(Interval(3,5) <= Interval(3,6));
  CHECK(Interval(3,6) <= Interval(3,6));
  CHECK_FALSE(Interval(3,8) <= Interval(3,6));
  CHECK_FALSE(Interval(4,5) <= Interval(3,6));
  CHECK_FALSE(Interval(4,6) <= Interval(3,6));
  CHECK_FALSE(Interval(4,8) <= Interval(3,6));
  CHECK_FALSE(Interval(6,8) <= Interval(3,6));
  CHECK_FALSE(Interval(7,8) <= Interval(3,6));
}

TEST_CASE("interval::operator>") {
  CHECK_FALSE(Interval(1,2) > Interval(3,6));
  CHECK_FALSE(Interval(1,3) > Interval(3,6));
  CHECK_FALSE(Interval(1,4) > Interval(3,6));
  CHECK_FALSE(Interval(1,6) > Interval(3,6));
  CHECK_FALSE(Interval(1,8) > Interval(3,6));
  CHECK_FALSE(Interval(3,5) > Interval(3,6));
  CHECK_FALSE(Interval(3,6) > Interval(3,6));
  CHECK(Interval(3,8) > Interval(3,6));
  CHECK(Interval(4,5) > Interval(3,6));
  CHECK(Interval(4,6) > Interval(3,6));
  CHECK(Interval(4,8) > Interval(3,6));
  CHECK(Interval(6,8) > Interval(3,6));
  CHECK(Interval(7,8) > Interval(3,6));
}

TEST_CASE("interval::operator>=") {
  CHECK_FALSE(Interval(1,2) >= Interval(3,6));
  CHECK_FALSE(Interval(1,3) >= Interval(3,6));
  CHECK_FALSE(Interval(1,4) >= Interval(3,6));
  CHECK_FALSE(Interval(1,6) >= Interval(3,6));
  CHECK_FALSE(Interval(1,8) >= Interval(3,6));
  CHECK_FALSE(Interval(3,5) >= Interval(3,6));
  CHECK(Interval(3,6) >= Interval(3,6));
  CHECK(Interval(3,8) >= Interval(3,6));
  CHECK(Interval(4,5) >= Interval(3,6));
  CHECK(Interval(4,6) >= Interval(3,6));
  CHECK(Interval(4,8) >= Interval(3,6));
  CHECK(Interval(6,8) >= Interval(3,6));
  CHECK(Interval(7,8) >= Interval(3,6));
}

TEST_CASE("interval::sgn") {
  CHECK(sgn(Interval(-3,-2)) == -1);
  CHECK(sgn(Interval(-3,0)) == -1);
  CHECK(sgn(Interval(-3,2)) == 0);
  CHECK(sgn(Interval(0,2)) == 1);
  CHECK(sgn(Interval(1,2)) == 1);
}

TEST_CASE("interval::pow") {
  CHECK(pow(Interval(1,2), 2) == Interval(1,4));
  CHECK(pow(Interval(1,2), 3) == Interval(1,8));
  CHECK(pow(Interval(-1,1), 5) == Interval(-1,1));
  CHECK(pow(Interval(-1,1), 4) == Interval(0, false, 1, true));
}

TEST_CASE("interval::operator+") {
  CHECK(Interval(1,2) + Interval(1,2) == Interval(2,4));
  CHECK(Interval(1,2) + Interval(-1,2) == Interval(0,4));
  CHECK(Interval(-1,1) + Interval(1,2) == Interval(0,3));
  CHECK(Interval(-1,1) + Interval(-1,2) == Interval(-2,3));
}

TEST_CASE("interval::operator*") {
  CHECK(Interval(1,2) * Interval(1,2) == Interval(1,4));
  CHECK(Interval(1,2) * Interval(-1,2) == Interval(-2,4));
  CHECK(Interval(-1,1) * Interval(1,2) == Interval(-2,2));
  CHECK(Interval(-1,1) * Interval(-1,2) == Interval(-2,2));
}

TEST_CASE("interval::operator<<") {
    Interval i(1, 3);
    std::stringstream out;
    out << i;
    CHECK(out.str() == "(1, 3)");
}
