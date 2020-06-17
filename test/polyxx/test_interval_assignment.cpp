#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("interval_assignment") {
  IntervalAssignment a;

  Variable x("x");
  Variable y("y");
  Variable z("z");

  CHECK_FALSE(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(x, Interval(1, 2));

  CHECK(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(y, Interval(2, 3));

  CHECK(a.has(x));
  CHECK(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(z, Interval(3, 4));

  CHECK(a.has(x));
  CHECK(a.has(y));
  CHECK(a.has(z));

  CHECK(a.get(x) == Interval(1, 2));
  CHECK(a.get(y) == Interval(2, 3));
  CHECK(a.get(z) == Interval(3, 4));

  a.unset(x);
  a.unset(z);

  CHECK_FALSE(a.has(x));
  CHECK(a.has(y));
  CHECK_FALSE(a.has(z));

  a.clear();

  CHECK_FALSE(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(x, Interval(4, 5));
  a.set(z, Interval(5, 6));

  CHECK(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK(a.has(z));
  CHECK(a.get(x) == Interval(4, 5));
  CHECK(a.get(z) == Interval(5, 6));
}