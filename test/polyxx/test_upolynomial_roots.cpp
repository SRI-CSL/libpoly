#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("upolynomial::count_real_roots") {
  {
    UPolynomial p = UPolynomial({-2, 0, 1})*UPolynomial({-2, 0, 1}) * UPolynomial({-3, 0, 1});
    CHECK(count_real_roots(p, RationalInterval(-5,5)) == 4);
    CHECK(count_real_roots(p, RationalInterval(-5,0)) == 2);
  }
}

TEST_CASE("upolynomial::isolate_real_roots") {
  {
    UPolynomial p = UPolynomial({-2, 0, 1}) * UPolynomial({-3, 0, 1});
    std::vector<AlgebraicNumber> roots = isolate_real_roots(p);
    CHECK(roots.size() == 4);
    CHECK(roots[0] == AlgebraicNumber(UPolynomial({-3, 0, 1}), DyadicInterval(-2, -1)));
    CHECK(roots[1] == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1)));
    CHECK(roots[2] == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2)));
    CHECK(roots[3] == AlgebraicNumber(UPolynomial({-3, 0, 1}), DyadicInterval(1, 2)));
  }
  {
    UPolynomial p = UPolynomial({-2, 0, 1})*UPolynomial({-2, 0, 1}) * UPolynomial({-3, 0, 1});
    std::vector<AlgebraicNumber> roots = isolate_real_roots(p);
    CHECK(roots.size() == 4);
    CHECK(roots[0] == AlgebraicNumber(UPolynomial({-3, 0, 1}), DyadicInterval(-2, -1)));
    CHECK(roots[1] == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(-2, -1)));
    CHECK(roots[2] == AlgebraicNumber(UPolynomial({-2, 0, 1}), DyadicInterval(1, 2)));
    CHECK(roots[3] == AlgebraicNumber(UPolynomial({-3, 0, 1}), DyadicInterval(1, 2)));
  }
}
