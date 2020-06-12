#pragma once

#include "../integer.h"
#include "utils.h"

#include <iosfwd>

namespace poly {

  class Integer;

  /** Represents an integer ring. Use IntegerRing::Z for the ring of all
   * integers. */
  class IntegerRing {
    /** The actual ring. */
    deleting_unique_ptr<lp_int_ring_t> mRing;

    /** Construct the ring of all integers Z. */
    IntegerRing();

   public:
    /** Construct the ring with the given modulus. */
    IntegerRing(const Integer& m, bool is_prime);

    /** Get a non-const pointer to the internal lp_int_ring_t. Handle with care!
     */
    lp_int_ring_t* get_internal();
    /** Get a const pointer to the internal lp_int_ring_t. */
    const lp_int_ring_t* get_internal() const;

    /** The ring of all integers Z. */
    static IntegerRing Z;
  };

  /** Compare two rings for equality. */
  bool operator==(const IntegerRing& lhs, const IntegerRing& rhs);
  /** Compare two rings for inequality. */
  bool operator!=(const IntegerRing& lhs, const IntegerRing& rhs);

  /** Stream the given IntegerRing to an output stream. */
  std::ostream& operator<<(std::ostream& os, const IntegerRing& ir);

}  // namespace poly
