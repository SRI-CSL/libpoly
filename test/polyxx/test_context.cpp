#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("context::discriminant") {
  Context ctx;
  Variable x(ctx, "x");
  Polynomial p1(ctx, x);
  Polynomial d = discriminant(p1);
  Polynomial p2 = p1 + d;
  CHECK(1 == d);
  CHECK(p2 == Polynomial(ctx, x) + 1);
}

TEST_CASE("context::integer") {
  Context ctx;
  Polynomial p(ctx);
  Polynomial p1 = p + 1;
  Polynomial p2 = p + Integer(1);
  CHECK(p1 == p2);
}
