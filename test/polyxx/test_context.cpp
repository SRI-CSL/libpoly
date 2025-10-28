#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("context::discriminant") {
  Context ctx;
  Variable x(ctx, "x");
  Polynomial p1(ctx, x);
  Polynomial d = poly::discriminant(p1); // d built in the default context
  Polynomial p2 = p1 + d; // assertion failure!
}
