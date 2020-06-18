#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

using namespace poly;

TEST_CASE("value::constructors") {
    CHECK(is_none(Value()));
    CHECK(is_integer(Value(13)));
    CHECK(is_algebraic_number(Value(AlgebraicNumber(DyadicRational(1)))));
    CHECK(is_dyadic_rational(Value(DyadicRational(1))));
    CHECK(is_integer(Value(Integer(13))));
    CHECK(is_rational(Value(Rational(27,4))));
    CHECK(is_plus_infinity(Value(Value::plus_infty())));
    CHECK(is_minus_infinity(Value(Value::minus_infty())));
}

TEST_CASE("value::operator==") {
    CHECK(Value() == Value());

    CHECK(Value::minus_infty() == Value::minus_infty());
    CHECK_FALSE(Value::minus_infty() == Value(-20));
    CHECK_FALSE(Value::minus_infty() == Value(Integer(-15)));
    CHECK_FALSE(Value::minus_infty() == Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value::minus_infty() == Value(Rational(1,2)));
    CHECK_FALSE(Value::minus_infty() == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value::minus_infty() == Value::plus_infty());

    CHECK(Value(-20) == Value(-20));
    CHECK_FALSE(Value(-20) == Value(Integer(-15)));
    CHECK_FALSE(Value(-20) == Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(-20) == Value(Rational(1,2)));
    CHECK_FALSE(Value(-20) == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(-20) == Value::plus_infty());

    CHECK(Value(Integer(-15)) == Value(Integer(-15)));
    CHECK_FALSE(Value(Integer(-15)) == Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(Integer(-15)) == Value(Rational(1,2)));
    CHECK_FALSE(Value(Integer(-15)) == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Integer(-15)) == Value::plus_infty());

    CHECK(Value(DyadicRational(-6,2)) == Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) == Value(Rational(1,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) == Value::plus_infty());

    CHECK(Value(Rational(1,2)) == Value(Rational(1,2)));
    CHECK_FALSE(Value(Rational(1,2)) == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Rational(1,2)) == Value::plus_infty());

    CHECK(Value(AlgebraicNumber(5)) == Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(AlgebraicNumber(5)) == Value::plus_infty());

    CHECK(Value::plus_infty() == Value::plus_infty());
}

TEST_CASE("value::operator!=") {
    CHECK_FALSE(Value() != Value());
    
    CHECK_FALSE(Value::minus_infty() != Value::minus_infty());
    CHECK(Value::minus_infty() != Value(-20));
    CHECK(Value::minus_infty() != Value(Integer(-15)));
    CHECK(Value::minus_infty() != Value(DyadicRational(-6,2)));
    CHECK(Value::minus_infty() != Value(Rational(1,2)));
    CHECK(Value::minus_infty() != Value(AlgebraicNumber(5)));
    CHECK(Value::minus_infty() != Value::plus_infty());

    CHECK_FALSE(Value(-20) != Value(-20));
    CHECK(Value(-20) != Value(Integer(-15)));
    CHECK(Value(-20) != Value(DyadicRational(-6,2)));
    CHECK(Value(-20) != Value(Rational(1,2)));
    CHECK(Value(-20) != Value(AlgebraicNumber(5)));
    CHECK(Value(-20) != Value::plus_infty());

    CHECK_FALSE(Value(Integer(-15)) != Value(Integer(-15)));
    CHECK(Value(Integer(-15)) != Value(DyadicRational(-6,2)));
    CHECK(Value(Integer(-15)) != Value(Rational(1,2)));
    CHECK(Value(Integer(-15)) != Value(AlgebraicNumber(5)));
    CHECK(Value(Integer(-15)) != Value::plus_infty());

    CHECK_FALSE(Value(DyadicRational(-6,2)) != Value(DyadicRational(-6,2)));
    CHECK(Value(DyadicRational(-6,2)) != Value(Rational(1,2)));
    CHECK(Value(DyadicRational(-6,2)) != Value(AlgebraicNumber(5)));
    CHECK(Value(DyadicRational(-6,2)) != Value::plus_infty());

    CHECK_FALSE(Value(Rational(1,2)) != Value(Rational(1,2)));
    CHECK(Value(Rational(1,2)) != Value(AlgebraicNumber(5)));
    CHECK(Value(Rational(1,2)) != Value::plus_infty());

    CHECK_FALSE(Value(AlgebraicNumber(5)) != Value(AlgebraicNumber(5)));
    CHECK(Value(AlgebraicNumber(5)) != Value::plus_infty());

    CHECK_FALSE(Value::plus_infty() != Value::plus_infty());
}

TEST_CASE("value::operator<") {
    CHECK_FALSE(Value() < Value());
    
    CHECK_FALSE(Value::minus_infty() < Value::minus_infty());
    CHECK(Value::minus_infty() < Value(-20));
    CHECK(Value::minus_infty() < Value(Integer(-15)));
    CHECK(Value::minus_infty() < Value(DyadicRational(-6,2)));
    CHECK(Value::minus_infty() < Value(Rational(1,2)));
    CHECK(Value::minus_infty() < Value(AlgebraicNumber(5)));
    CHECK(Value::minus_infty() < Value::plus_infty());

    CHECK_FALSE(Value(-20) < Value(-20));
    CHECK(Value(-20) < Value(Integer(-15)));
    CHECK(Value(-20) < Value(DyadicRational(-6,2)));
    CHECK(Value(-20) < Value(Rational(1,2)));
    CHECK(Value(-20) < Value(AlgebraicNumber(5)));
    CHECK(Value(-20) < Value::plus_infty());

    CHECK_FALSE(Value(Integer(-15)) < Value(Integer(-15)));
    CHECK(Value(Integer(-15)) < Value(DyadicRational(-6,2)));
    CHECK(Value(Integer(-15)) < Value(Rational(1,2)));
    CHECK(Value(Integer(-15)) < Value(AlgebraicNumber(5)));
    CHECK(Value(Integer(-15)) < Value::plus_infty());

    CHECK_FALSE(Value(DyadicRational(-6,2)) < Value(DyadicRational(-6,2)));
    CHECK(Value(DyadicRational(-6,2)) < Value(Rational(1,2)));
    CHECK(Value(DyadicRational(-6,2)) < Value(AlgebraicNumber(5)));
    CHECK(Value(DyadicRational(-6,2)) < Value::plus_infty());

    CHECK_FALSE(Value(Rational(1,2)) < Value(Rational(1,2)));
    CHECK(Value(Rational(1,2)) < Value(AlgebraicNumber(5)));
    CHECK(Value(Rational(1,2)) < Value::plus_infty());

    CHECK_FALSE(Value(AlgebraicNumber(5)) < Value(AlgebraicNumber(5)));
    CHECK(Value(AlgebraicNumber(5)) < Value::plus_infty());

    CHECK_FALSE(Value::plus_infty() < Value::plus_infty());
}

TEST_CASE("value::operator<=") {
    CHECK(Value() <= Value());
    
    CHECK(Value::minus_infty() <= Value::minus_infty());
    CHECK(Value::minus_infty() <= Value(-20));
    CHECK(Value::minus_infty() <= Value(Integer(-15)));
    CHECK(Value::minus_infty() <= Value(DyadicRational(-6,2)));
    CHECK(Value::minus_infty() <= Value(Rational(1,2)));
    CHECK(Value::minus_infty() <= Value(AlgebraicNumber(5)));
    CHECK(Value::minus_infty() <= Value::plus_infty());

    CHECK(Value(-20) <= Value(-20));
    CHECK(Value(-20) <= Value(Integer(-15)));
    CHECK(Value(-20) <= Value(DyadicRational(-6,2)));
    CHECK(Value(-20) <= Value(Rational(1,2)));
    CHECK(Value(-20) <= Value(AlgebraicNumber(5)));
    CHECK(Value(-20) <= Value::plus_infty());

    CHECK(Value(Integer(-15)) <= Value(Integer(-15)));
    CHECK(Value(Integer(-15)) <= Value(DyadicRational(-6,2)));
    CHECK(Value(Integer(-15)) <= Value(Rational(1,2)));
    CHECK(Value(Integer(-15)) <= Value(AlgebraicNumber(5)));
    CHECK(Value(Integer(-15)) <= Value::plus_infty());

    CHECK(Value(DyadicRational(-6,2)) <= Value(DyadicRational(-6,2)));
    CHECK(Value(DyadicRational(-6,2)) <= Value(Rational(1,2)));
    CHECK(Value(DyadicRational(-6,2)) <= Value(AlgebraicNumber(5)));
    CHECK(Value(DyadicRational(-6,2)) <= Value::plus_infty());

    CHECK(Value(Rational(1,2)) <= Value(Rational(1,2)));
    CHECK(Value(Rational(1,2)) <= Value(AlgebraicNumber(5)));
    CHECK(Value(Rational(1,2)) <= Value::plus_infty());

    CHECK(Value(AlgebraicNumber(5)) <= Value(AlgebraicNumber(5)));
    CHECK(Value(AlgebraicNumber(5)) <= Value::plus_infty());

    CHECK(Value::plus_infty() <= Value::plus_infty());
}

TEST_CASE("value::operator>") {
    CHECK_FALSE(Value() > Value());
    
    CHECK_FALSE(Value::minus_infty() > Value::minus_infty());
    CHECK_FALSE(Value::minus_infty() > Value(-20));
    CHECK_FALSE(Value::minus_infty() > Value(Integer(-15)));
    CHECK_FALSE(Value::minus_infty() > Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value::minus_infty() > Value(Rational(1,2)));
    CHECK_FALSE(Value::minus_infty() > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value::minus_infty() > Value::plus_infty());

    CHECK_FALSE(Value(-20) > Value(-20));
    CHECK_FALSE(Value(-20) > Value(Integer(-15)));
    CHECK_FALSE(Value(-20) > Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(-20) > Value(Rational(1,2)));
    CHECK_FALSE(Value(-20) > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(-20) > Value::plus_infty());

    CHECK_FALSE(Value(Integer(-15)) > Value(Integer(-15)));
    CHECK_FALSE(Value(Integer(-15)) > Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(Integer(-15)) > Value(Rational(1,2)));
    CHECK_FALSE(Value(Integer(-15)) > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Integer(-15)) > Value::plus_infty());

    CHECK_FALSE(Value(DyadicRational(-6,2)) > Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) > Value(Rational(1,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) > Value::plus_infty());

    CHECK_FALSE(Value(Rational(1,2)) > Value(Rational(1,2)));
    CHECK_FALSE(Value(Rational(1,2)) > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Rational(1,2)) > Value::plus_infty());

    CHECK_FALSE(Value(AlgebraicNumber(5)) > Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(AlgebraicNumber(5)) > Value::plus_infty());

    CHECK_FALSE(Value::plus_infty() > Value::plus_infty());
}

TEST_CASE("value::operator>=") {
    CHECK(Value() >= Value());
    
    CHECK(Value::minus_infty() >= Value::minus_infty());
    CHECK_FALSE(Value::minus_infty() >= Value(-20));
    CHECK_FALSE(Value::minus_infty() >= Value(Integer(-15)));
    CHECK_FALSE(Value::minus_infty() >= Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value::minus_infty() >= Value(Rational(1,2)));
    CHECK_FALSE(Value::minus_infty() >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value::minus_infty() >= Value::plus_infty());

    CHECK(Value(-20) >= Value(-20));
    CHECK_FALSE(Value(-20) >= Value(Integer(-15)));
    CHECK_FALSE(Value(-20) >= Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(-20) >= Value(Rational(1,2)));
    CHECK_FALSE(Value(-20) >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(-20) >= Value::plus_infty());


    CHECK(Value(Integer(-15)) >= Value(Integer(-15)));
    CHECK_FALSE(Value(Integer(-15)) >= Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(Integer(-15)) >= Value(Rational(1,2)));
    CHECK_FALSE(Value(Integer(-15)) >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Integer(-15)) >= Value::plus_infty());

    CHECK(Value(DyadicRational(-6,2)) >= Value(DyadicRational(-6,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) >= Value(Rational(1,2)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(DyadicRational(-6,2)) >= Value::plus_infty());

    CHECK(Value(Rational(1,2)) >= Value(Rational(1,2)));
    CHECK_FALSE(Value(Rational(1,2)) >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(Rational(1,2)) >= Value::plus_infty());

    CHECK(Value(AlgebraicNumber(5)) >= Value(AlgebraicNumber(5)));
    CHECK_FALSE(Value(AlgebraicNumber(5)) >= Value::plus_infty());

    CHECK(Value::plus_infty() >= Value::plus_infty());
}

TEST_CASE("value::sgn") {
    CHECK(sgn(Value::minus_infty()) == -1);
    CHECK(sgn(Value(-20)) == -1);
    CHECK(sgn(Value(Integer(-15))) == -1);
    CHECK(sgn(Value(DyadicRational(-6,2))) == -1);
    CHECK(sgn(Value(Rational(1,2))) == 1);
    CHECK(sgn(Value(AlgebraicNumber(5))) == 1);
    CHECK(sgn(Value::plus_infty()) == 1);
}

TEST_CASE("value::value_between") {
    Value a(AlgebraicNumber(3));
    Value b(AlgebraicNumber(7));
    Value s = value_between(a, true, b, true);
    CHECK(a < s);
    CHECK(s < b);
}
