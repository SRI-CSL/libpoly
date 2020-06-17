#include "polyxx/integer_ring.h"

#include <iostream>
#include <utility>

#include "polyxx/integer.h"

namespace poly {

  IntegerRing IntegerRing::Z;

  /** A deleter for an std::unique_ptr holding a lp_int_ring_t pointer */
  void int_ring_deleter(lp_int_ring_t* ptr) { lp_int_ring_detach(ptr); }

  IntegerRing::IntegerRing() : mRing(lp_Z, int_ring_deleter) {}

  IntegerRing::IntegerRing(const Integer& m, bool is_prime)
      : mRing(lp_int_ring_create(m.get_internal(), is_prime ? 1 : 0),
              int_ring_deleter) {}

  lp_int_ring_t* IntegerRing::get_internal() { return mRing.get(); }
  const lp_int_ring_t* IntegerRing::get_internal() const { return mRing.get(); }

  bool operator==(const IntegerRing& lhs, const IntegerRing& rhs) {
    return lp_int_ring_equal(lhs.get_internal(), rhs.get_internal());
  }
  bool operator!=(const IntegerRing& lhs, const IntegerRing& rhs) {
    return !(lhs == rhs);
  }

  std::ostream& operator<<(std::ostream& os, const IntegerRing& ir) {
    return stream_ptr(os, lp_int_ring_to_string(ir.get_internal()));
  }

}  // namespace poly
