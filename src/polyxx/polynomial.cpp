#include "polyxx/polynomial.h"

#include "feasibility_set.h"
#include "poly.h"
#include "variable_list.h"

#include <cassert>
#include <utility>

namespace poly {

  namespace detail {
    /** Return the context of a polynomial. */
    const lp_polynomial_context_t* context(const Polynomial& p) {
      return lp_polynomial_get_context(p.get_internal());
    }
    /** Return the context of two polynomials. Asserts that the context objects
     * are equal. */
    const lp_polynomial_context_t* context(const Polynomial& lhs,
                                           const Polynomial& rhs) {
      assert(lp_polynomial_context_equal(context(lhs), context(rhs)));
      (void)rhs;
      return context(lhs);
    }
  }  // namespace detail

  /** A deleter for an std::unique_ptr holding a lp_polynomial_t pointer */
  void polynomial_deleter(lp_polynomial_t* ptr) { lp_polynomial_delete(ptr); }

  Polynomial::Polynomial(lp_polynomial_t* poly)
      : mPoly(poly, polynomial_deleter) {}
  Polynomial::Polynomial(const lp_polynomial_t* poly)
      : mPoly(lp_polynomial_new_copy(poly), polynomial_deleter) {}
  Polynomial::Polynomial(const lp_polynomial_context_t* c)
      : mPoly(lp_polynomial_new(c), polynomial_deleter) {}
  Polynomial::Polynomial(const Context& c)
      : Polynomial(c.get_polynomial_context()) {}
  Polynomial::Polynomial() : Polynomial(Context::get_context()) {}

  Polynomial::Polynomial(const lp_polynomial_context_t* c, const Variable& v)
      : mPoly(lp_polynomial_alloc(), polynomial_deleter) {
    assert(lp_variable_db_is_valid(c->var_db, v.get_internal()));
    lp_polynomial_construct_simple(get_internal(), c,
                                   Integer(1).get_internal(),
                                   v.get_internal(), 1);
  }
  Polynomial::Polynomial(const Context& c, const Variable& v)
      : Polynomial(c.get_polynomial_context(), v) {}
  Polynomial::Polynomial(const Variable& v) : Polynomial(Context::get_context(), v) {}

  Polynomial::Polynomial(const lp_polynomial_context_t* c, const Integer& i, const Variable& v, unsigned n)
      : mPoly(lp_polynomial_alloc(), polynomial_deleter) {
    assert(lp_variable_db_is_valid(c->var_db, v.get_internal()));
    lp_polynomial_construct_simple(get_internal(), c,
                                   i.get_internal(), v.get_internal(), n);
  }
  Polynomial::Polynomial(const Context& c, const Integer& i, const Variable& v, unsigned n)
      : Polynomial(c.get_polynomial_context(), i, v, n) {}
  Polynomial::Polynomial(const Integer& i, const Variable& v, unsigned n)
      : Polynomial(Context::get_context(), i, v, n) {}

  Polynomial::Polynomial(const lp_polynomial_context_t* c, const Integer & i)
      : mPoly(lp_polynomial_alloc(), polynomial_deleter) {
    lp_polynomial_construct_simple(get_internal(), c,
                                   i.get_internal(), lp_variable_null, 0);
  }
  Polynomial::Polynomial(const Context& c, const Integer& i) : Polynomial(c.get_polynomial_context(), i) {}
  Polynomial::Polynomial(const Integer& i) : Polynomial(Context::get_context(), i){};

  Polynomial::Polynomial(const Polynomial& p)
      : mPoly(lp_polynomial_new_copy(p.get_internal()), polynomial_deleter) {}
  Polynomial::Polynomial(Polynomial&& p) noexcept
      : mPoly(std::move(p.mPoly)) {}

  Polynomial& Polynomial::operator=(const Polynomial& p) {
    mPoly.reset(lp_polynomial_new_copy(p.get_internal()));
    return *this;
  }
  Polynomial& Polynomial::operator=(Polynomial&& p) noexcept {
    mPoly = std::move(p.mPoly);
    return *this;
  }

  lp_polynomial_t* Polynomial::get_internal() { return mPoly.get(); }
  const lp_polynomial_t* Polynomial::get_internal() const {
    return mPoly.get();
  }
  lp_polynomial_t* Polynomial::release() { return mPoly.release(); }

  void swap(Polynomial& lhs, Polynomial& rhs) {
    lp_polynomial_swap(lhs.get_internal(), rhs.get_internal());
  }

  std::size_t hash(const Polynomial& p) {
    return lp_polynomial_hash(p.get_internal());
  }

  std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
    return stream_ptr(os, lp_polynomial_to_string(p.get_internal()));
  }

  bool is_zero(const Polynomial& p) {
    return lp_polynomial_is_zero(p.get_internal());
  }
  bool is_constant(const Polynomial& p) {
    return lp_polynomial_is_constant(p.get_internal());
  }
  bool is_linear(const Polynomial& p) {
    return lp_polynomial_is_linear(p.get_internal());
  }
  bool is_lc_constant(const Polynomial& p) {
    return lp_polynomial_lc_is_constant(p.get_internal());
  }
  int lc_sgn(const Polynomial& p) {
    return lp_polynomial_lc_sgn(p.get_internal());
  }
  std::size_t degree(const Polynomial& p) {
    return lp_polynomial_degree(p.get_internal());
  }
  Variable main_variable(const Polynomial& p) {
    return Variable(lp_polynomial_top_variable(p.get_internal()));
  }
  Polynomial coefficient(const Polynomial& p, std::size_t k) {
    Polynomial res(detail::context(p));
    lp_polynomial_get_coefficient(res.get_internal(), p.get_internal(), k);
    return res;
  }
  Polynomial leading_coefficient(const Polynomial& p) {
    return coefficient(p, degree(p));
  }
  std::vector<Polynomial> coefficients(const Polynomial& p) {
    std::vector<Polynomial> res;
    for (std::size_t deg = 0; deg <= degree(p); ++deg) {
      auto coeff = coefficient(p, deg);
      if (lp_polynomial_is_constant(coeff.get_internal())) continue;
      res.emplace_back(coeff);
    }
    return res;
  }

  bool is_univariate(const Polynomial& p) {
    return lp_polynomial_is_univariate(p.get_internal());
  }
  UPolynomial to_univariate(const Polynomial& p) {
    assert(is_univariate(p));
    return UPolynomial(lp_polynomial_to_univariate(p.get_internal()));
  }
  bool is_univariate_over_assignment(const Polynomial& p, const Assignment& a) {
    return lp_polynomial_is_univariate_m(p.get_internal(), a.get_internal());
  }
  bool is_assigned_over_assignment(const Polynomial& p, const Assignment& a) {
    return lp_polynomial_is_assigned(p.get_internal(), a.get_internal());
  }
  UPolynomial to_univariate(const Polynomial& p, const Assignment& a) {
    return UPolynomial(lp_polynomial_to_univariate_m(p.get_internal(), a.get_internal()));
  }
  int sgn(const Polynomial& p, const Assignment& a) {
    return lp_polynomial_sgn(p.get_internal(), a.get_internal());
  }
  Value evaluate(const Polynomial& p, const Assignment& a) {
    return Value(lp_polynomial_evaluate(p.get_internal(), a.get_internal()));
  }
  bool evaluate_constraint(const Polynomial& p, const Assignment& a,
                           SignCondition sc) {
    return lp_polynomial_constraint_evaluate(
        p.get_internal(), to_sign_condition(sc), a.get_internal());
  }
  int evaluate_constraint_subs(const Polynomial& p, const Assignment& a,
                               SignCondition sc) {
    return lp_polynomial_constraint_evaluate_subs(
      p.get_internal(), to_sign_condition(sc), a.get_internal());
  }
  Interval evaluate(const Polynomial& p, const IntervalAssignment& a) {
    Interval res;
    lp_polynomial_interval_value(p.get_internal(), a.get_internal(),
                                 res.get_internal());
    return res;
  }

  bool operator==(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_eq(lhs.get_internal(), rhs.get_internal());
  }
  bool operator==(const Polynomial& lhs, const Integer& rhs) {
    Polynomial tmp(detail::context(lhs), rhs);
    return lp_polynomial_eq(lhs.get_internal(), tmp.get_internal());
  }
  bool operator==(const Integer& lhs, const Polynomial& rhs) {
    Polynomial tmp(detail::context(rhs), lhs);
    return lp_polynomial_eq(tmp.get_internal(), rhs.get_internal());
  }
  bool operator!=(const Polynomial& lhs, const Polynomial& rhs) {
    return !(lhs == rhs);
  }
  bool operator!=(const Polynomial& lhs, const Integer& rhs) {
    return !(lhs == rhs);
  }
  bool operator!=(const Integer& lhs, const Polynomial& rhs) {
    return !(lhs == rhs);
  }
  bool operator<(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<(const Polynomial& lhs, const Integer& rhs) {
    Polynomial tmp(detail::context(lhs), rhs);
    return lp_polynomial_cmp(lhs.get_internal(), tmp.get_internal()) < 0;
  }
  bool operator<(const Integer& lhs, const Polynomial& rhs) {
    Polynomial tmp(detail::context(rhs), lhs);
    return lp_polynomial_cmp(tmp.get_internal(), rhs.get_internal()) < 0;
  }
  bool operator<=(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator<=(const Polynomial& lhs, const Integer& rhs) {
    Polynomial tmp(detail::context(lhs), rhs);
    return lp_polynomial_cmp(lhs.get_internal(), tmp.get_internal()) <= 0;
  }
  bool operator<=(const Integer& lhs, const Polynomial& rhs) {
    Polynomial tmp(detail::context(rhs), lhs);
    return lp_polynomial_cmp(tmp.get_internal(), rhs.get_internal()) <= 0;
  }
  bool operator>(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>(const Polynomial& lhs, const Integer& rhs) {
    Polynomial tmp(detail::context(lhs), rhs);
    return lp_polynomial_cmp(lhs.get_internal(), tmp.get_internal()) > 0;
  }
  bool operator>(const Integer& lhs, const Polynomial& rhs) {
    Polynomial tmp(detail::context(rhs), lhs);
    return lp_polynomial_cmp(tmp.get_internal(), rhs.get_internal()) > 0;
  }
  bool operator>=(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
  }
  bool operator>=(const Polynomial& lhs, const Integer& rhs) {
    Polynomial tmp(detail::context(lhs), rhs);
    return lp_polynomial_cmp(lhs.get_internal(), tmp.get_internal()) >= 0;
  }
  bool operator>=(const Integer& lhs, const Polynomial& rhs) {
    Polynomial tmp(detail::context(rhs), lhs);
    return lp_polynomial_cmp(tmp.get_internal(), rhs.get_internal()) >= 0;
  }

  Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_add(res.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return res;
  }
  Polynomial operator+(const Polynomial& lhs, const Integer& rhs) {
    Polynomial res(lhs);
    res += rhs;
    return res;
  }
  Polynomial operator+(const Integer& lhs, const Polynomial& rhs) {
    Polynomial res(rhs);
    res += lhs;
    return res;
  }
  Polynomial& operator+=(Polynomial& lhs, const Polynomial& rhs) {
    lp_polynomial_add(lhs.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return lhs;
  }
  Polynomial& operator+=(Polynomial& lhs, const Integer& rhs) {
    lp_monomial_t monomial;
    lp_monomial_construct(detail::context(lhs), &monomial);
    lp_monomial_set_coefficient(detail::context(lhs), &monomial,
                                rhs.get_internal());
    lp_polynomial_add_monomial(lhs.get_internal(), &monomial);
    lp_monomial_destruct(&monomial);
    return lhs;
  }

  Polynomial& add_mul(Polynomial& lhs, const Polynomial& rhs1,
                      const Polynomial& rhs2) {
    lp_polynomial_add_mul(lhs.get_internal(), rhs1.get_internal(),
                          rhs2.get_internal());
    return lhs;
  }

  Polynomial operator-(const Polynomial& p) {
    Polynomial res(detail::context(p));
    lp_polynomial_neg(res.get_internal(), p.get_internal());
    return res;
  }
  Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_sub(res.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return res;
  }
  Polynomial operator-(const Polynomial& lhs, const Integer& rhs) {
    return lhs + (-rhs);
  }
  Polynomial operator-(const Integer& lhs, const Polynomial& rhs) {
    return (-rhs) + lhs;
  }
  Polynomial& operator-=(Polynomial& lhs, const Polynomial& rhs) {
    lp_polynomial_sub(lhs.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return lhs;
  }
  Polynomial& operator-=(Polynomial& lhs, const Integer& rhs) {
    lhs += (-rhs);
    return lhs;
  }

  Polynomial& sub_mul(Polynomial& lhs, const Polynomial& rhs1,
                      const Polynomial& rhs2) {
    lp_polynomial_sub_mul(lhs.get_internal(), rhs1.get_internal(),
                          rhs2.get_internal());
    return lhs;
  }

  Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_mul(res.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return res;
  }
  Polynomial operator*(const Polynomial& lhs, const Integer& rhs) {
    Polynomial res(detail::context(lhs));
    lp_polynomial_mul_integer(res.get_internal(), lhs.get_internal(),
                              rhs.get_internal());
    return res;
  }
  Polynomial operator*(const Integer& lhs, const Polynomial& rhs) {
    return rhs * lhs;
  }
  Polynomial& operator*=(Polynomial& lhs, const Polynomial& rhs) {
    lp_polynomial_mul(lhs.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return lhs;
  }
  Polynomial& operator*=(Polynomial& lhs, const Integer& rhs) {
    lp_polynomial_mul_integer(lhs.get_internal(), lhs.get_internal(),
                              rhs.get_internal());
    return lhs;
  }

  Polynomial shl(const Polynomial& lhs, unsigned exp) {
    Polynomial res(detail::context(lhs));
    lp_polynomial_shl(res.get_internal(), lhs.get_internal(), exp);
    return res;
  }
  Polynomial pow(const Polynomial& lhs, unsigned exp) {
    Polynomial res(detail::context(lhs));
    lp_polynomial_pow(res.get_internal(), lhs.get_internal(), exp);
    return res;
  }

  bool divides(const Polynomial& lhs, const Polynomial& rhs) {
    return lp_polynomial_divides(lhs.get_internal(), rhs.get_internal());
  }
  Polynomial div(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_div(res.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return res;
  }
  Polynomial rem(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_rem(res.get_internal(), lhs.get_internal(),
                      rhs.get_internal());
    return res;
  }
  Polynomial prem(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_prem(res.get_internal(), lhs.get_internal(),
                       rhs.get_internal());
    return res;
  }
  Polynomial sprem(const Polynomial& lhs, const Polynomial& rhs) {
    Polynomial res(detail::context(lhs, rhs));
    lp_polynomial_sprem(res.get_internal(), lhs.get_internal(),
                        rhs.get_internal());
    return res;
  }
  std::pair<Polynomial, Polynomial> p_div_rem(const Polynomial& lhs,
                                              const Polynomial& rhs) {
    Polynomial d(detail::context(lhs, rhs));
    Polynomial r(detail::context(lhs, rhs));
    lp_polynomial_pdivrem(d.get_internal(), r.get_internal(), lhs.get_internal(),
                         rhs.get_internal());
    return {d, r};
  }
  std::pair<Polynomial, Polynomial> sp_div_rem(const Polynomial& lhs,
                                               const Polynomial& rhs) {
    Polynomial d(detail::context(lhs, rhs));
    Polynomial r(detail::context(lhs, rhs));
    lp_polynomial_spdivrem(d.get_internal(), r.get_internal(), lhs.get_internal(),
                          rhs.get_internal());
    return {d, r};
  }
  std::pair<Polynomial, Polynomial> div_rem(const Polynomial& lhs,
                                            const Polynomial& rhs) {
    Polynomial d(detail::context(lhs, rhs));
    Polynomial r(detail::context(lhs, rhs));
    lp_polynomial_divrem(d.get_internal(), r.get_internal(), lhs.get_internal(),
                         rhs.get_internal());
    return {d, r};
  }

  Polynomial derivative(const Polynomial& p) {
    Polynomial res(detail::context(p));
    lp_polynomial_derivative(res.get_internal(), p.get_internal());
    return res;
  }

  Polynomial gcd(const Polynomial& p, const Polynomial& q) {
    Polynomial res(detail::context(p, q));
    lp_polynomial_gcd(res.get_internal(), p.get_internal(), q.get_internal());
    return res;
  }
  Polynomial lcm(const Polynomial& p, const Polynomial& q) {
    Polynomial res(detail::context(p, q));
    lp_polynomial_lcm(res.get_internal(), p.get_internal(), q.get_internal());
    return res;
  }

  Polynomial content(const Polynomial& p) {
    Polynomial res(detail::context(p));
    lp_polynomial_cont(res.get_internal(), p.get_internal());
    return res;
  }
  Polynomial primitive_part(const Polynomial& p) {
    Polynomial res(detail::context(p));
    lp_polynomial_pp(res.get_internal(), p.get_internal());
    return res;
  }
  std::pair<Polynomial, Polynomial> content_primitive_part(
      const Polynomial& p) {
    Polynomial cn(detail::context(p));
    Polynomial pp(detail::context(p));
    lp_polynomial_pp_cont(pp.get_internal(), cn.get_internal(),
                          p.get_internal());
    return {cn, pp};
  }

  Polynomial resultant(const Polynomial& p, const Polynomial& q) {
    Polynomial res(detail::context(p, q));
    lp_polynomial_resultant(res.get_internal(), p.get_internal(),
                            q.get_internal());
    return res;
  }
  Polynomial discriminant(const Polynomial& p) {
    if (degree(p) == 1) {
      // Derivative is constant, making the resultant trivial
      // (and resultant() does not cope with that)
      return Polynomial(detail::context(p), Integer(1));
    }
    return div(resultant(p, derivative(p)), leading_coefficient(p));
  }

  std::vector<Polynomial> psc(const Polynomial& p, const Polynomial& q) {
    std::size_t size = std::min(degree(p), degree(q)) + 1;
    std::vector<lp_polynomial_t*> tmp;
    for (std::size_t i = 0; i < size; ++i) {
      tmp.push_back(lp_polynomial_new(detail::context(p, q)));
    }
    lp_polynomial_psc(tmp.data(), p.get_internal(), q.get_internal());
    std::vector<Polynomial> res;
    for (std::size_t i = 0; i < size; ++i) {
      res.emplace_back(tmp[i]);
    }
    tmp.clear();
    return res;
  }

  std::vector<Polynomial> subres(const Polynomial& p, const Polynomial& q) {
    std::size_t size = std::min(degree(p), degree(q)) + 1;
    std::vector<lp_polynomial_t*> tmp;
    for (std::size_t i = 0; i < size; ++i) {
      tmp.push_back(lp_polynomial_new(detail::context(p, q)));
    }
    lp_polynomial_subres(tmp.data(), p.get_internal(), q.get_internal());
    std::vector<Polynomial> res;
    for (std::size_t i = 0; i < size; ++i) {
      res.emplace_back(tmp[i]);
    }
    tmp.clear();
    return res;
  }

  std::vector<Polynomial> square_free_factors(const Polynomial& p) {
    lp_polynomial_t** factors = nullptr;
    std::size_t* multiplicities = nullptr;
    std::size_t size = 0;
    lp_polynomial_factor_square_free(p.get_internal(), &factors,
                                     &multiplicities, &size);

    std::vector<Polynomial> res;
    for (std::size_t i = 0; i < size; ++i) {
      res.emplace_back(factors[i]);
    }
    free(factors);
    free(multiplicities);

    return res;
  }

  std::vector<Polynomial> content_free_factors(const Polynomial& p) {
    lp_polynomial_t** factors = nullptr;
    std::size_t* multiplicities = nullptr;
    std::size_t size = 0;
    lp_polynomial_factor_content_free(p.get_internal(), &factors,
                                      &multiplicities, &size);

    std::vector<Polynomial> res;
    for (std::size_t i = 0; i < size; ++i) {
      res.emplace_back(factors[i]);
    }
    free(factors);
    free(multiplicities);

    return res;
  }

  std::vector<Value> isolate_real_roots(const Polynomial& p,
                                        const Assignment& a) {
    lp_value_t* roots = new lp_value_t[degree(p)];
    std::size_t roots_size;
    lp_polynomial_roots_isolate(p.get_internal(), a.get_internal(), roots,
                                &roots_size);
    std::vector<Value> res;
    for (std::size_t i = 0; i < roots_size; ++i) {
      res.emplace_back();
      lp_value_construct_copy(res.back().get_internal(), &roots[i]);
    }
    for (std::size_t i = 0; i < roots_size; ++i) {
      lp_value_destruct(&roots[i]);
    }
    delete[] roots;
    return res;
  }

  std::vector<Interval> infeasible_regions(const Polynomial& p,
                                           const Assignment& a,
                                           SignCondition sc) {
    lp_feasibility_set_t* feasible = lp_polynomial_constraint_get_feasible_set(
        p.get_internal(), to_sign_condition(sc), 0, a.get_internal());

    std::vector<Interval> regions;

    Value last_value = Value::minus_infty();
    int last_open = 0;

    for (std::size_t i = 0; i < feasible->size; ++i) {
      const lp_interval_t& cur = feasible->intervals[i];
      Value lower(&cur.a);

      if (lower.get_internal()->type == LP_VALUE_MINUS_INFINITY) {
        // Do nothing if we start at -infty.
      } else if (last_value < lower) {
        // There is an infeasible open interval
        regions.emplace_back(last_value, !last_open, lower, !cur.a_open);
      } else if (last_open && cur.a_open && last_value == lower) {
        // There is an infeasible point interval
        regions.emplace_back(last_value);
      }
      if (cur.is_point) {
        last_value = std::move(lower);
        last_open = false;
      } else {
        last_value = Value(&cur.b);
        last_open = cur.b_open;
      }
    }

    if (last_value.get_internal()->type != LP_VALUE_PLUS_INFINITY) {
      // Add missing interval to +infty
      regions.emplace_back(last_value, !last_open, Value::plus_infty(), true);
    }

    lp_feasibility_set_delete(feasible);

    return regions;
  }

}  // namespace poly
