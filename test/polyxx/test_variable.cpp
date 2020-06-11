#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("variable") {
  Variable v1;
  Variable v2("foo");

  CHECK(v1 == v1);
  CHECK_FALSE(v1 == v2);
  CHECK(v2 == v2);
  CHECK_FALSE(v1 != v1);
  CHECK(v1 != v2);
  CHECK_FALSE(v2 != v2);
}