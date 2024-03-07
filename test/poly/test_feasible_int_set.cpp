#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../polyxx/doctest.h"

#include <poly.h>
#include <feasibility_set_int.h>

#include <polyxx.h>
using namespace poly;
using namespace std;

static
vector<lp_integer_t>* integer_list_gen(std::initializer_list<int> values) {
  auto result = new vector<lp_integer_t>();
  for (int x : values) {
    lp_integer_t tmp;
    lp_integer_construct_from_int(lp_Z, &tmp, x);
    result->push_back(tmp);
  }
  return result;
}

static
void integer_list_del(vector<lp_integer_t> *list) {
  for (lp_integer_t i : *list) {
    lp_integer_destruct(&i);
  }
  delete list;
}

TEST_CASE("feasibility_set_int::construct") {
  lp_integer_t prime;
  lp_integer_construct_from_int(lp_Z, &prime, 13);
  lp_int_ring_t *K = lp_int_ring_create(&prime, 1);

  SUBCASE("feasibility_set_int::construct::full") {
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_full(K);
    CHECK(lp_feasibility_set_int_is_full(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::empty") {
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_empty(K);
    CHECK(lp_feasibility_set_int_is_empty(set));
    lp_feasibility_set_int_delete(set);
  }
  SUBCASE("feasibility_set_int::construct::integer") {
    auto *ints = integer_list_gen({2,1,4,3});
    lp_feasibility_set_int_t *set = lp_feasibility_set_int_new_from_integer(K, ints->data(), ints->size(), false);
    CHECK_FALSE(lp_feasibility_set_int_is_empty(set));
    CHECK_FALSE(lp_feasibility_set_int_is_full(set));
    CHECK(lp_feasibility_set_int_size_approx(set) == 4);
    lp_feasibility_set_int_delete(set);
    integer_list_del(ints);
  }

  CHECK_EQ(K->ref_count, 1);
  lp_int_ring_detach(K);
  lp_integer_destruct(&prime);
}
