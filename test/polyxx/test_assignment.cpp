#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("assignment") {
  Assignment a;

  Variable x("x");
  Variable y("y");
  Variable z("z");

  CHECK_FALSE(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(x, Value(1));

  CHECK(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(y, Value(2));

  CHECK(a.has(x));
  CHECK(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(z, Value(3));

  CHECK(a.has(x));
  CHECK(a.has(y));
  CHECK(a.has(z));

  CHECK(a.get(x) == Value(1));
  CHECK(a.get(y) == Value(2));
  CHECK(a.get(z) == Value(3));

  a.unset(x);
  a.unset(z);

  CHECK_FALSE(a.has(x));
  CHECK(a.has(y));
  CHECK_FALSE(a.has(z));

  a.clear();

  CHECK_FALSE(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK_FALSE(a.has(z));

  a.set(x, Value(4));
  a.set(z, Value(5));

  CHECK(a.has(x));
  CHECK_FALSE(a.has(y));
  CHECK(a.has(z));
  CHECK(a.get(x) == Value(4));
  CHECK(a.get(z) == Value(5));

  std::stringstream out;
  out << a;
  CHECK(out.str() == "[x -> 4, z -> 5]");
}
