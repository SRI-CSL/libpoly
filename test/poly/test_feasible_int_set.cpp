#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../polyxx/doctest.h"

#include <poly.h>
#include <feasibility_set_int.h>

#include <polyxx.h>
using namespace poly;
using namespace std;

/**
 * Struct to handle a list of lp_integer_t (and proper deallocate it)
 */
struct integer_list {
private:
  vector<lp_integer_t> *v;
public:
  integer_list() : v(new vector<lp_integer_t>()) {}
  integer_list(initializer_list<int> values) : integer_list() {
    for (int x : values) {
      lp_integer_t tmp;
      lp_integer_construct_from_int(lp_Z, &tmp, x);
      this->v->push_back(tmp);
    }
  }
  integer_list(const integer_list &other) : integer_list() {
    for (const lp_integer_t x : *other.v) {
      lp_integer_t tmp;
      lp_integer_construct_copy(lp_Z, &tmp, &x);
      this->v->push_back(tmp);
    }
  }
  integer_list(integer_list &&other) noexcept : v(other.v) {
    other.v = nullptr;
  }
  ~integer_list() {
    if (this->v == nullptr) return;
    for (lp_integer_t i : *this->v) {
      lp_integer_destruct(&i);
    }
    delete this->v;
  }
  lp_feasibility_set_int_t* gen_set(lp_int_ring_t *K, bool inverted) const {
    return lp_feasibility_set_int_new_from_integer(K, this->v->data(), this->v->size(), inverted);
  }
  lp_feasibility_set_int_t* gen_set(IntegerRing &K, bool inverted) const {
    return this->gen_set(K.get_internal(), inverted);
  }
  const lp_integer_t* data() const { return this->v->data(); }
  size_t size() const { return this->v->size(); }
};


TEST_CASE("feasibility_set_int::construct") {
  lp_integer_t prime;
  lp_integer_construct_from_int(lp_Z, &prime, 13);
  lp_int_ring_t *K = lp_int_ring_create(&prime, 1);

  SUBCASE("feasibility_set_int::construct::full") {
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_full(K);
    CHECK(!lp_feasibility_set_int_is_empty(set));
    CHECK(lp_feasibility_set_int_is_full(set));
    CHECK(!lp_feasibility_set_int_is_point(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::empty") {
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_empty(K);
    CHECK(lp_feasibility_set_int_is_empty(set));
    CHECK(!lp_feasibility_set_int_is_full(set));
    CHECK(!lp_feasibility_set_int_is_point(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::integer") {
    integer_list l({2,1,4,3});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, l.data(), l.size(), false);
    CHECK_FALSE(lp_feasibility_set_int_is_empty(set));
    CHECK_FALSE(lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_size_approx(set) == 4);
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::point") {
    integer_list l({2});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, l.data(), l.size(), false);
    CHECK(!lp_feasibility_set_int_is_empty(set));
    CHECK(!lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_is_point(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::point_inv") {
    integer_list l({0,1,/*2,*/3,4,5,6,7,8,9,10,11,12});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, l.data(), l.size(), true);
    CHECK(!lp_feasibility_set_int_is_empty(set));
    CHECK(!lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_is_point(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::integer::invert") {
    integer_list l({2,1,4,3});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, l.data(), l.size(), true);
    CHECK_FALSE(lp_feasibility_set_int_is_empty(set));
    CHECK_FALSE(lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_size_approx(set) == 9);
    CHECK(!lp_feasibility_set_int_contains(set, Integer(3).get_internal()));
    CHECK(!lp_feasibility_set_int_contains(set, Integer(1).get_internal()));
    CHECK(lp_feasibility_set_int_contains(set, Integer(-1).get_internal()));
    CHECK(lp_feasibility_set_int_contains(set, Integer(13).get_internal()));
    CHECK(!lp_feasibility_set_int_contains(set, Integer(14).get_internal()));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::integer::duplicates") {
    integer_list l({1,2,2,14});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, l.data(), l.size(), false);
    CHECK_FALSE(lp_feasibility_set_int_is_empty(set));
    CHECK_FALSE(lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_size_approx(set) == 2);
    lp_feasibility_set_int_delete(set);
  }

  CHECK_EQ(K->ref_count, 1);
  lp_int_ring_detach(K);
  lp_integer_destruct(&prime);
}

TEST_CASE("feasibility_set_int::intersect::inv") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set2 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::intersect::inv1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }

  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::intersect::inv_dual") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set2 = integer_list({1,2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::intersect::inv_dual1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv_dual2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv_dual") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::inv4_dual") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }

  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::union::inv") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set2 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::union::inv1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::union::inv_dual") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set2 = integer_list({1,2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::union::inv_dual1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv_dual2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv_dual3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union::inv_dual4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set2));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::union") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set2 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::union1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 3);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::intersect") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set2 = integer_list({2,3}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *set2_inv = integer_list({0,1,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::intersect1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 1);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 1);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 1);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, set2_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 1);
    CHECK(status == LP_FEASIBILITY_SET_INT_NEW);
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(set2);
  lp_feasibility_set_int_delete(set2_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::intersect_empty") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *empty = lp_feasibility_set_int_new_empty(K.get_internal());

  SUBCASE("feasibility_set_int::intersect_empty1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, empty, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 0);
    CHECK(status == LP_FEASIBILITY_SET_INT_EMPTY);
    CHECK(lp_feasibility_set_int_is_empty(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect_empty2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(empty, set1, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 0);
    CHECK(status == LP_FEASIBILITY_SET_INT_EMPTY);
    CHECK(lp_feasibility_set_int_is_empty(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect_empty3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, empty, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 0);
    CHECK(status == LP_FEASIBILITY_SET_INT_EMPTY);
    CHECK(lp_feasibility_set_int_is_empty(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect_empty4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(empty, set1_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 0);
    CHECK(status == LP_FEASIBILITY_SET_INT_EMPTY);
    CHECK(lp_feasibility_set_int_is_empty(s));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(empty);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::intersect_full") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *full = lp_feasibility_set_int_new_full(K.get_internal());

  SUBCASE("feasibility_set_int::intersect::intersect_full1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1, full, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::intersect_full2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(full, set1, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::intersect_full3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(set1_inv, full, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::intersect::intersect_full4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_intersect_with_status(full, set1_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(full);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::union_full") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *full = lp_feasibility_set_int_new_full(K.get_internal());

  SUBCASE("feasibility_set_int::union_full1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, full, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 7);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_is_full(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_full2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(full, set1, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 7);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_is_full(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_full3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, full, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 7);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_is_full(s));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_full4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(full, set1_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 7);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_is_full(s));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(full);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::union_empty") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *empty = lp_feasibility_set_int_new_empty(K.get_internal());

  SUBCASE("feasibility_set_int::union_empty1") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1, empty, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_empty2") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(empty, set1, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_empty3") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(set1_inv, empty, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S1);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  SUBCASE("feasibility_set_int::union_empty4") {
    lp_feasibility_set_int_status_t status;
    lp_feasibility_set_int_t *s = lp_feasibility_set_int_union_with_status(empty, set1_inv, &status);
    CHECK(lp_feasibility_set_int_size_approx(s) == 2);
    CHECK(status == LP_FEASIBILITY_SET_INT_S2);
    CHECK(lp_feasibility_set_int_eq(s, set1));
    lp_feasibility_set_int_delete(s);
  }
  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(empty);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::size") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);
  lp_feasibility_set_int_t *empty = lp_feasibility_set_int_new_empty(K.get_internal());
  lp_feasibility_set_int_t *full = lp_feasibility_set_int_new_full(K.get_internal());

  SUBCASE("feasibility_set_int::size::1") {
    Integer size;
    lp_feasibility_set_int_size(set1, size.get_internal());
    size_t size_approx = lp_feasibility_set_int_size_approx(set1);
    CHECK(size == 2);
    CHECK(size_approx == 2);
  }
  SUBCASE("feasibility_set_int::size::2") {
    Integer size;
    lp_feasibility_set_int_size(set1_inv, size.get_internal());
    size_t size_approx = lp_feasibility_set_int_size_approx(set1_inv);
    CHECK(size == 2);
    CHECK(size_approx == 2);
  }
  SUBCASE("feasibility_set_int::size::full") {
    Integer size;
    lp_feasibility_set_int_size(full, size.get_internal());
    size_t size_approx = lp_feasibility_set_int_size_approx(full);
    CHECK(size == prime);
    CHECK(size_approx == prime);
  }
  SUBCASE("feasibility_set_int::size::empty") {
    Integer size;
    lp_feasibility_set_int_size(empty, size.get_internal());
    size_t size_approx = lp_feasibility_set_int_size_approx(empty);
    CHECK(size == 0);
    CHECK(size_approx == 0);
  }

  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);
  lp_feasibility_set_int_delete(full);
  lp_feasibility_set_int_delete(empty);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}

TEST_CASE("feasibility_set_int::pick") {
  Integer prime(7);
  IntegerRing K(prime, true);

  lp_feasibility_set_int_t *set1 = integer_list({1,2}).gen_set(K, false);
  lp_feasibility_set_int_t *set1_inv = integer_list({0,3,4,5,6}).gen_set(K, true);

  SUBCASE("feasibility_set_int::pick") {
    Integer pick;
    lp_feasibility_set_int_pick_value(set1, pick.get_internal());

    CHECK((pick == 1 || pick == 2));
  }
  SUBCASE("feasibility_set_int::pick_inv") {
    Integer pick;
    lp_feasibility_set_int_pick_value(set1_inv, pick.get_internal());
    CHECK((pick == 1 || pick == 2));
  }

  lp_feasibility_set_int_delete(set1);
  lp_feasibility_set_int_delete(set1_inv);

  CHECK_EQ(K.get_internal()->ref_count, 1);
}
