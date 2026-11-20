#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("polynomial::properties") {
  Variable x("x");
  Variable y("y");
  Polynomial p1;
  Polynomial p2(Integer(5));
  Polynomial p3 = 3*x;
  Polynomial p4 = 3*x*y*y + 7*x;

  CHECK(is_zero(p1));
  CHECK_FALSE(is_zero(p2));
  CHECK_FALSE(is_zero(p3));
  CHECK_FALSE(is_zero(p4));

  CHECK(is_constant(p1));
  CHECK(is_constant(p2));
  CHECK_FALSE(is_constant(p3));
  CHECK_FALSE(is_constant(p4));

  CHECK_FALSE(is_linear(p1));
  CHECK_FALSE(is_linear(p2));
  CHECK(is_linear(p3));
  CHECK_FALSE(is_linear(p4));

  CHECK(is_lc_constant(p1));
  CHECK(is_lc_constant(p2));
  CHECK(is_lc_constant(p3));
  CHECK_FALSE(is_lc_constant(p4));

  CHECK(degree(p1) == 0);
  CHECK(degree(p2) == 0);
  CHECK(degree(p3) == 1);
  CHECK(degree(p4) == 2);

  CHECK(is_univariate(p3));
}

TEST_CASE("polynomial::resultant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  Polynomial r = resultant(p, q);
  CHECK(r == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
}

TEST_CASE("polynomial::discriminant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial d = discriminant(p);
  CHECK(d == 11337408 * pow(y, 5) - 35095680 * pow(y, 4) + 34197120 * pow(y, 3) - 14999040 * pow(y, 2) + 3099840 * y - 246656);
}
TEST_CASE("polynomial::psc") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  auto tmp = psc(p, q);

  CHECK(tmp.size() == 6);
  CHECK(tmp[0] == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
  CHECK(tmp[1] == 416745 * y * y * y - 416745 * y * y + 138915 * y - 15435);
  CHECK(tmp[2] == Integer(0));
  CHECK(tmp[3] == Integer(0));
  CHECK(tmp[4] == Integer(-45));
  CHECK(tmp[5] == Integer(7));
}

TEST_CASE("polynomial::operator<<") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  std::stringstream out;
  out << p;
  CHECK(out.str() == "1*x^6 + 2*x^5 + (3*y - 1)");
}

TEST_CASE("polynomial::constructor") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  lp_polynomial_t* ptr = p.get_internal();

  Polynomial p_cp(p);
  Polynomial p_mv(std::move(p));
  CHECK(ptr == p_mv.get_internal());
  CHECK(ptr != p_cp.get_internal());
  CHECK(p_cp == p_mv);

  Polynomial p_mv_a, p_cp_a;
  p_cp_a = p_cp;
  p_mv_a = std::move(p_mv);
  CHECK(p_cp == p_cp_a);
  CHECK(ptr == p_mv_a.get_internal());
}

TEST_CASE("polynomial::evaluate") {
  Variable x("x");
  Variable y("y");
  Polynomial p1 = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial p2 = 1 * y * pow(x + 1, 4);

  // full assignment
  Assignment a;
  a.set(x, -1);
  a.set(y, 2);
  CHECK(a.has(x));
  CHECK(a.has(y));

  CHECK(evaluate(p1, a) == Value(4));
  CHECK(evaluate_constraint(p1, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p1, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p1, a, SignCondition::EQ) == false);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::EQ) == 0);

  CHECK(evaluate(p2, a) == Value(0));
  CHECK(evaluate_constraint(p2, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p2, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p2, a, SignCondition::EQ) == true);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::EQ) == 1);

  // partial assignment
  a.unset(y);
  CHECK(a.has(x));
  CHECK(!a.has(y));

  CHECK(evaluate_constraint_subs(p1, a, SignCondition::GE) == -1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::LT) == -1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::EQ) == -1);

  CHECK(evaluate_constraint(p2, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p2, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p2, a, SignCondition::EQ) == true);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::EQ) == 1);
}

TEST_CASE("polynomial::infeasible_regions") {
  Variable x("x");
  Variable y("y");
  Polynomial p = pow(y, 2) - x;
  Assignment a;

  // y^2 - 4 has the roots -2 and 2
  a.set(x, 4);
  std::vector<Value> roots;
  std::vector<Interval> regions = infeasible_regions(p, a, SignCondition::LT, roots);
  REQUIRE(roots.size() == 2);
  CHECK(roots[0] == Value(-2));
  CHECK(roots[1] == Value(2));
  REQUIRE(regions.size() == 2);
  CHECK(regions[0] == Interval(Value::minus_infty(), true, -2, false));
  CHECK(regions[1] == Interval(2, false, Value::plus_infty(), true));
  CHECK(regions == infeasible_regions(p, a, SignCondition::LT));

  // the roots are appended, the regions are point intervals at the roots
  regions = infeasible_regions(p, a, SignCondition::NE, roots);
  REQUIRE(roots.size() == 4);
  CHECK(roots[0] == Value(-2));
  CHECK(roots[1] == Value(2));
  CHECK(roots[2] == Value(-2));
  CHECK(roots[3] == Value(2));
  REQUIRE(regions.size() == 2);
  CHECK(regions[0] == Interval(Value(-2)));
  CHECK(regions[1] == Interval(Value(2)));
  CHECK(regions == infeasible_regions(p, a, SignCondition::NE));

  // y^2 - 2 has the roots -sqrt(2) and sqrt(2)
  a.set(x, 2);
  roots.clear();
  regions = infeasible_regions(p, a, SignCondition::GE, roots);
  CHECK(roots == isolate_real_roots(p, a));
  REQUIRE(roots.size() == 2);
  CHECK(is_algebraic_number(roots[0]));
  CHECK(is_algebraic_number(roots[1]));
  CHECK(roots[0] < roots[1]);
  REQUIRE(regions.size() == 1);
  CHECK(regions[0] == Interval(roots[0], true, roots[1], true));
  CHECK(regions == infeasible_regions(p, a, SignCondition::GE));

  // y^2 has a double root at 0
  a.set(x, 0);
  roots.clear();
  regions = infeasible_regions(p, a, SignCondition::LE, roots);
  REQUIRE(roots.size() == 1);
  CHECK(roots[0] == Value(0));
  REQUIRE(regions.size() == 2);
  CHECK(regions[0] == Interval(Value::minus_infty(), true, 0, true));
  CHECK(regions[1] == Interval(0, true, Value::plus_infty(), true));
  CHECK(regions == infeasible_regions(p, a, SignCondition::LE));

  // y^2 + 1 has no real roots
  a.set(x, -1);
  roots.clear();
  regions = infeasible_regions(p, a, SignCondition::LT, roots);
  CHECK(roots.empty());
  REQUIRE(regions.size() == 1);
  CHECK(is_full(regions[0]));
  CHECK(infeasible_regions(p, a, SignCondition::GT, roots).empty());
  CHECK(roots.empty());
}

TEST_CASE("polynomial::infeasible_regions constant under assignment") {
  Variable x("x");
  Variable y("y");
  // evaluates to the constant 3 for x = 1
  Polynomial p = (x - 1) * y + 3;
  Assignment a;
  a.set(x, 1);

  std::vector<Value> roots;
  std::vector<Interval> regions = infeasible_regions(p, a, SignCondition::LT, roots);
  CHECK(roots.empty());
  REQUIRE(regions.size() == 1);
  CHECK(is_full(regions[0]));

  CHECK(infeasible_regions(p, a, SignCondition::GT, roots).empty());
  CHECK(roots.empty());
}
