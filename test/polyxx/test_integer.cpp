#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <polyxx.h>

using namespace poly;

TEST_CASE("integer::constructors") {
    Integer i;
    CHECK(i == Integer(0l));
}

TEST_CASE("integer::bit_size") {
    CHECK(bit_size(Integer(5)) == 3);
    CHECK(bit_size(Integer(64)) == 7);
}

TEST_CASE("integer::comparisons") {
    CHECK(Integer(1) == Integer(1));
    CHECK_FALSE(Integer(1) == Integer(2));
    CHECK_FALSE(Integer(2) == Integer(1));

    CHECK_FALSE(Integer(1) != Integer(1));
    CHECK(Integer(1) != Integer(2));
    CHECK(Integer(2) != Integer(1));

    CHECK_FALSE(Integer(1) < Integer(1));
    CHECK(Integer(1) < Integer(2));
    CHECK_FALSE(Integer(2) < Integer(1));

    CHECK(Integer(1) <= Integer(1));
    CHECK(Integer(1) <= Integer(2));
    CHECK_FALSE(Integer(2) <= Integer(1));

    CHECK_FALSE(Integer(1) > Integer(1));
    CHECK_FALSE(Integer(1) > Integer(2));
    CHECK(Integer(2) > Integer(1));

    CHECK(Integer(1) >= Integer(1));
    CHECK_FALSE(Integer(1) >= Integer(2));
    CHECK(Integer(2) >= Integer(1));
}

TEST_CASE("integer::divides") {
    CHECK_FALSE(divides(Integer(15), Integer(5)));
    CHECK(divides(Integer(5), Integer(15)));
}

TEST_CASE("integer::swap") {
    Integer a(5);
    Integer b(10);
    swap(a, b);
    CHECK(a == Integer(10));
    CHECK(b == Integer(5));
}

TEST_CASE("integer::++") {
    Integer a(5);
    CHECK(++a == Integer(6));
    CHECK(a == Integer(6));
    CHECK(a++ == Integer(6));
    CHECK(a == Integer(7));
}

TEST_CASE("integer::--") {
    Integer a(5);
    CHECK(--a == Integer(4));
    CHECK(a == Integer(4));
    CHECK(a-- == Integer(4));
    CHECK(a == Integer(3));
}

TEST_CASE("integer::operator+") {
    CHECK(Integer(3) + Integer(4) == Integer(7));
    Integer a(3);
    a += Integer(4);
    CHECK(a == Integer(7));
}

TEST_CASE("integer::operator-") {
    CHECK(Integer(3) - Integer(2) == Integer(1));
    Integer a(3);
    a -= Integer(2);
    CHECK(a == Integer(1));
    CHECK(-Integer(3) == Integer(-3));
}

TEST_CASE("integer::abs") {
    CHECK(abs(Integer(-1)) == Integer(1));
    CHECK(abs(Integer()) == Integer());
    CHECK(abs(Integer(1)) == Integer(1));
}

TEST_CASE("integer::operator*") {
    CHECK(Integer(3) * Integer(4) == Integer(12));
    CHECK(Integer(3) * 4 == Integer(12));
    CHECK(3 * Integer(4) == Integer(12));
    Integer a(3);
    a *= Integer(4);
    CHECK(a == Integer(12));
    a *= 4;
    CHECK(a == Integer(48));
}

TEST_CASE("integer::mul_pow2") {
    CHECK(mul_pow2(Integer(3), 0) == Integer(3));
    CHECK(mul_pow2(Integer(3), 1) == Integer(6));
    CHECK(mul_pow2(Integer(3), 2) == Integer(12));
}

TEST_CASE("integer::pow") {
    CHECK(pow(Integer(3), 0) == Integer(1));
    CHECK(pow(Integer(3), 1) == Integer(3));
    CHECK(pow(Integer(3), 2) == Integer(9));
}

TEST_CASE("integer::sqrt") {
    CHECK(sqrt(Integer(0)) == Integer(0));
    CHECK(sqrt(Integer(1)) == Integer(1));
    CHECK(sqrt(Integer(2)) == Integer(1));
    CHECK(sqrt(Integer(3)) == Integer(1));
    CHECK(sqrt(Integer(4)) == Integer(2));
    CHECK(sqrt(Integer(15)) == Integer(3));
    CHECK(sqrt(Integer(16)) == Integer(4));
    CHECK(sqrt(Integer(17)) == Integer(4));
}

TEST_CASE("integer::add_mul") {
    Integer a(3);
    add_mul(a, Integer(2), Integer(3));
    CHECK(a == Integer(9));
    add_mul(a, Integer(2), 3);
    CHECK(a == Integer(15));
}

TEST_CASE("integer::sub_mul") {
    Integer a(20);
    sub_mul(a, Integer(2), Integer(3));
    CHECK(a == Integer(14));
}

TEST_CASE("integer::operator/") {
    CHECK(Integer(13) / Integer(4) == Integer(3));
    CHECK(Integer(12) / Integer(4) == Integer(3));
    CHECK(Integer(11) / Integer(4) == Integer(2));
}

TEST_CASE("integer::operator/=") {
    Integer a(20);
    a /= Integer(5);
    CHECK(a == Integer(4));
    a /= Integer(3);
    CHECK(a == Integer(1));
}

TEST_CASE("integer::operator%") {
    CHECK(Integer(13) % Integer(4) == Integer(1));
    CHECK(Integer(12) % Integer(4) == Integer(0));
    CHECK(Integer(11) % Integer(4) == Integer(3));
}

TEST_CASE("integer::operator%=") {
    Integer a(23);
    a %= Integer(5);
    CHECK(a == Integer(3));
    a %= Integer(3);
    CHECK(a == Integer(0));
}

TEST_CASE("integer::div_exact") {
    CHECK(div_exact(Integer(16), Integer(4)) == Integer(4));
    CHECK(div_exact(Integer(12), Integer(4)) == Integer(3));
    CHECK(div_exact(Integer(8), Integer(4)) == Integer(2));
}

TEST_CASE("integer::div_exact") {
    Integer rem;
    CHECK(div_rem(rem, Integer(13), Integer(4)) == Integer(3));
    CHECK(rem == Integer(1));
    CHECK(div_rem(rem, Integer(12), Integer(4)) == Integer(3));
    CHECK(rem == Integer(0));
    CHECK(div_rem(rem, Integer(11), Integer(4)) == Integer(2));
    CHECK(rem == Integer(3));
}

TEST_CASE("integer::div_exact") {
    Integer rem;
    CHECK(div_rem_pow2(rem, Integer(13), 2) == Integer(3));
    CHECK(rem == Integer(1));
    CHECK(div_rem_pow2(rem, Integer(12), 2) == Integer(3));
    CHECK(rem == Integer(0));
    CHECK(div_rem_pow2(rem, Integer(11), 2) == Integer(2));
    CHECK(rem == Integer(3));
}

TEST_CASE("integer::to_int") {
    CHECK(to_int(Integer(-10)) == -10);
    CHECK(to_int(Integer()) == 0);
    CHECK(to_int(Integer(5)) == 5);
}

TEST_CASE("integer::to_double") {
    CHECK(to_double(Integer(-10)) == -10.0);
    CHECK(to_double(Integer()) == 0.0);
    CHECK(to_double(Integer(5)) == 5.0);
}

TEST_CASE("integer::is_prime") {
    CHECK(is_prime(Integer(-3)));
    CHECK(is_prime(Integer(-2)));
    CHECK_FALSE(is_prime(Integer(-1)));
    CHECK_FALSE(is_prime(Integer()));
    CHECK_FALSE(is_prime(Integer(1)));
    CHECK(is_prime(Integer(2)));
    CHECK(is_prime(Integer(3)));
    CHECK_FALSE(is_prime(Integer(4)));
    CHECK(is_prime(Integer(5)));
    CHECK_FALSE(is_prime(Integer(6)));
    CHECK(is_prime(Integer(7)));
}

TEST_CASE("integer::is_zero") {
    CHECK_FALSE(is_zero(Integer(-1)));
    CHECK(is_zero(Integer()));
    CHECK_FALSE(is_zero(Integer(3)));
}

TEST_CASE("integer::sgn") {
    CHECK(sgn(Integer(-1)) == -1);
    CHECK(sgn(Integer()) == 0);
    CHECK(sgn(Integer(1)) == 1);
}

TEST_CASE("integer::gcd") {
    CHECK(gcd(Integer(15), Integer(125)) == Integer(5));
}

TEST_CASE("integer::lcm") {
    CHECK(lcm(Integer(15), Integer(35)) == Integer(105));
}

TEST_CASE("integer::operator<<") {
    Integer i(5);
    std::stringstream out;
    out << i;
    CHECK(out.str() == "5");
}
