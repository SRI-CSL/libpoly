#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("upolynomial::degree") {
  CHECK(degree(UPolynomial()) == 0);
  CHECK(degree(UPolynomial({1})) == 0);
  CHECK(degree(UPolynomial({1, 1})) == 1);
  CHECK(degree(UPolynomial({-2, 0, 1})) == 2);
}

TEST_CASE("upolynomial::leading_coefficient") {
  CHECK(leading_coefficient(UPolynomial()) == Integer(0));
  CHECK(leading_coefficient(UPolynomial({1})) == Integer(1));
  CHECK(leading_coefficient(UPolynomial({1, 2})) == Integer(2));
  CHECK(leading_coefficient(UPolynomial({-2, 0, 3})) == Integer(3));
}

TEST_CASE("upolynomial::constant_coefficient") {
  CHECK(constant_coefficient(UPolynomial()) == Integer(0));
  CHECK(constant_coefficient(UPolynomial({1})) == Integer(1));
  CHECK(constant_coefficient(UPolynomial({2, 2})) == Integer(2));
  CHECK(constant_coefficient(UPolynomial({2, 0, 3})) == Integer(2));
}

TEST_CASE("upolynomial::constant_coefficient") {
  CHECK(coefficients(UPolynomial()) == std::vector<Integer>({Integer(0)}));
  CHECK(coefficients(UPolynomial({1})) == std::vector<Integer>({Integer(1)}));
  CHECK(coefficients(UPolynomial({2, 2})) ==
        std::vector<Integer>({Integer(2), Integer(2)}));
  CHECK(coefficients(UPolynomial({2, 0, 3})) ==
        std::vector<Integer>({Integer(2), Integer(), Integer(3)}));
}

TEST_CASE("upolynomial::is_zero") {
  CHECK(is_zero(UPolynomial()));
  CHECK_FALSE(is_zero(UPolynomial({-1})));
  CHECK(is_zero(UPolynomial({0})));
  CHECK_FALSE(is_zero(UPolynomial({1})));
  CHECK_FALSE(is_zero(UPolynomial({2, 2})));
  CHECK_FALSE(is_zero(UPolynomial({0, 2})));
}

TEST_CASE("upolynomial::is_one") {
  CHECK_FALSE(is_one(UPolynomial()));
  CHECK_FALSE(is_one(UPolynomial({-1})));
  CHECK_FALSE(is_one(UPolynomial({0})));
  CHECK(is_one(UPolynomial({1})));
  CHECK_FALSE(is_one(UPolynomial({2, 2})));
  CHECK_FALSE(is_one(UPolynomial({0, 2})));
}

TEST_CASE("upolynomial::is_monic") {
  CHECK_FALSE(is_monic(UPolynomial()));
  CHECK_FALSE(is_monic(UPolynomial({-1})));
  CHECK_FALSE(is_monic(UPolynomial({0})));
  CHECK(is_monic(UPolynomial({1})));
  CHECK(is_monic(UPolynomial({0, 1})));
  CHECK_FALSE(is_monic(UPolynomial({0, 2})));
  CHECK(is_monic(UPolynomial({1, 1})));
  CHECK_FALSE(is_monic(UPolynomial({1, 2})));
  CHECK(is_monic(UPolynomial({2, 1})));
  CHECK_FALSE(is_monic(UPolynomial({2, 2})));
}

TEST_CASE("upolynomial::is_primitive") {
  CHECK_FALSE(is_primitive(UPolynomial({-1})));
  CHECK(is_primitive(UPolynomial({1})));
  CHECK(is_primitive(UPolynomial({0, 1})));
  CHECK_FALSE(is_primitive(UPolynomial({0, 2})));
  CHECK(is_primitive(UPolynomial({1, 1})));
  CHECK(is_primitive(UPolynomial({1, 2})));
  CHECK(is_primitive(UPolynomial({2, 1})));
  CHECK_FALSE(is_primitive(UPolynomial({2, 2})));
}

TEST_CASE("upolynomial::evaluate_at") {
  UPolynomial p({1, 1, 1});
  CHECK(evaluate_at(p, Integer(-1)) == Integer(1));
  CHECK(evaluate_at(p, Integer()) == Integer(1));
  CHECK(evaluate_at(p, Integer(1)) == Integer(3));
  CHECK(evaluate_at(p, Integer(2)) == Integer(7));

  CHECK(evaluate_at(p, Rational(-1, 2)) == Rational(3, 4));
  CHECK(evaluate_at(p, Rational()) == Rational(1));
  CHECK(evaluate_at(p, Rational(1, 3)) == Rational(13, 9));
  CHECK(evaluate_at(p, Rational(2, 3)) == Rational(19, 9));

  CHECK(evaluate_at(p, DyadicRational(-1, 1)) == DyadicRational(3, 2));
  CHECK(evaluate_at(p, DyadicRational()) == DyadicRational(1));
  CHECK(evaluate_at(p, DyadicRational(1, 2)) == DyadicRational(21, 4));
  CHECK(evaluate_at(p, DyadicRational(1, 3)) == DyadicRational(73, 6));
}

TEST_CASE("upolynomial::sign_at") {
  UPolynomial p({2, 1, -1});
  CHECK(sign_at(p, Integer(-2)) == -1);
  CHECK(sign_at(p, Integer(-1)) == 0);
  CHECK(sign_at(p, Integer()) == 1);
  CHECK(sign_at(p, Integer(1)) == 1);
  CHECK(sign_at(p, Integer(2)) == 0);

  CHECK(sign_at(p, Rational(-5, 2)) == -1);
  CHECK(sign_at(p, Rational(-1)) == 0);
  CHECK(sign_at(p, Rational(-1, 2)) == 1);
  CHECK(sign_at(p, Rational(2)) == 0);
  CHECK(sign_at(p, Rational(7, 3)) == -1);

  CHECK(sign_at(p, DyadicRational(-13, 2)) == -1);
  CHECK(sign_at(p, DyadicRational(-8, 3)) == 0);
  CHECK(sign_at(p, DyadicRational(-1, 5)) == 1);
  CHECK(sign_at(p, DyadicRational(8, 2)) == 0);
  CHECK(sign_at(p, DyadicRational(71, 3)) == -1);
}

TEST_CASE("upolynomial::subst_x_neg") {
  CHECK(subst_x_neg(UPolynomial({2, 1, -1})) == UPolynomial({2, -1, -1}));
}

TEST_CASE("upolynomial::operator-") {
  CHECK(-UPolynomial({2, 1, -1}) == UPolynomial({-2, -1, 1}));
}

TEST_CASE("upolynomial::neg") {
  UPolynomial p({2, 1, -1});
  neg(p);
  CHECK(p == UPolynomial({-2, -1, 1}));
}

TEST_CASE("upolynomial::operator+") {
  CHECK(UPolynomial({2, 1, -1}) + UPolynomial({1, 2, 3}) ==
        UPolynomial({3, 3, 2}));
}
TEST_CASE("upolynomial::operator-") {
  CHECK(UPolynomial({2, 1, -1}) - UPolynomial({1, 2, 3}) ==
        UPolynomial({1, -1, -4}));
}
TEST_CASE("upolynomial::operator*") {
  CHECK(UPolynomial({2, 1, -1}) * UPolynomial({1, 2, 3}) ==
        UPolynomial({2, 5, 7, 1, -3}));
  CHECK(UPolynomial({2, 1, -1}) * Integer(3) == UPolynomial({6, 3, -3}));
  CHECK(Integer(3) * UPolynomial({2, 1, -1}) == UPolynomial({6, 3, -3}));
}

TEST_CASE("upolynomial::pow") {
  CHECK(pow(UPolynomial({2, 1, -1}), 3) ==
        UPolynomial({8, 12, -6, -11, 3, 3, -1}));
}

TEST_CASE("upolynomial::derivative") {
  CHECK(derivative(UPolynomial({2, 5, 7, 1, -3})) ==
        UPolynomial({5, 14, 3, -12}));
}

TEST_CASE("upolynomial::divides") {
  CHECK(divides(UPolynomial({2, 1, -1}), UPolynomial({2, 5, 7, 1, -3})));
  CHECK(divides(UPolynomial({1, 2, 3}), UPolynomial({2, 5, 7, 1, -3})));
}

TEST_CASE("upolynomial::div_degrees") {
  CHECK(div_degrees(UPolynomial({1, 0, 0, 2, 0, 0, 3}), 3) ==
        UPolynomial({1, 2, 3}));
}

TEST_CASE("upolynomial::div_exact") {
  CHECK(div_exact(UPolynomial({2, 5, 7, 1, -3}), UPolynomial({2, 1, -1})) ==
        UPolynomial({1, 2, 3}));
  CHECK(div_exact(UPolynomial({2, 5, 7, 1, -3}), UPolynomial({1, 2, 3})) ==
        UPolynomial({2, 1, -1}));
}

TEST_CASE("upolynomial::div_exact") {
  CHECK(div_exact(UPolynomial({2, 4, 6, 8}), Integer(2)) ==
        UPolynomial({1, 2, 3, 4}));
}

TEST_CASE("upolynomial::rem_exact") {
  CHECK(rem_exact(UPolynomial({2, 5, 7, 1, -3}), UPolynomial({3, 0, -1, 1})) ==
        UPolynomial({8, 14, 5}));
}

TEST_CASE("upolynomial::div_rem_exact") {
  auto res =
      div_rem_exact(UPolynomial({2, 5, 7, 1, -3}), UPolynomial({3, 0, -1, 1}));
  CHECK(res.first == UPolynomial({-2, -3}));
  CHECK(res.second == UPolynomial({8, 14, 5}));
}

TEST_CASE("upolynomial::div_rem_pseudo") {
  auto res =
      div_rem_pseudo(UPolynomial({2, 5, 7, 1, -3}), UPolynomial({3, 0, -1, 1}));
  CHECK(res.first == UPolynomial({-2, -3}));
  CHECK(res.second == UPolynomial({8, 14, 5}));
}

TEST_CASE("upolynomial::content") {
  CHECK(content(UPolynomial({12, 18, 24})) == Integer(6));
}

TEST_CASE("upolynomial::make_primitive") {
  UPolynomial p({12, 18, 24});
  make_primitive(p);
  CHECK(p == UPolynomial({2, 3, 4}));
}

TEST_CASE("upolynomial::primitive_part") {
  CHECK(primitive_part(UPolynomial({12, 18, 24})) == UPolynomial({2, 3, 4}));
}

TEST_CASE("upolynomial::gcd") {
  UPolynomial p({1, 2, 3, 4, 5});
  UPolynomial q({3, 4, 1});
  UPolynomial r({-3, 7, -9});
  CHECK(gcd(p * q, q * r) == q);
  CHECK(gcd(p * q, p * r) == p);
  CHECK(gcd(p * r, r * r) == -r);
}

TEST_CASE("upolynomial::square_free_factors") {
  UPolynomial p({1, 2, 3, 4, 5});
  auto factors = square_free_factors(p, true);
  UPolynomial prod(Integer(1));
  for (const auto& f : factors) {
    prod = f * prod;
  }
  CHECK(prod == p);
}

TEST_CASE("upolynomial::sturm_sequence") {
  auto seq = sturm_sequence(UPolynomial({2, 5, 7, 1, -3}));
  CHECK(seq.size() == 5);
  CHECK(seq[0] == UPolynomial({-2, -5, -7, -1, 3}));
  CHECK(seq[1] == UPolynomial({-5, -14, -3, 12}));
  CHECK(seq[2] == UPolynomial({101, 194, 171}));
  CHECK(seq[3] == UPolynomial({-733, 341}));
  CHECK(seq[4] == UPolynomial({-1}));
}

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
