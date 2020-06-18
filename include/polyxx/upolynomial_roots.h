#pragma once

#include "algebraic_number.h"
#include "rational_interval.h"
#include "upolynomial.h"

namespace poly {

  /** Count the real roots of a polynomial within a rational interval.*/
  std::size_t count_real_roots(const UPolynomial& p,
                               const RationalInterval& ri);

  /** Isolate the real roots of a UPolynomial, returning them as algebraic
   * numbers. The roots are sorted in increasing numerical order.
   */
  std::vector<AlgebraicNumber> isolate_real_roots(const UPolynomial& p);

}  // namespace poly
