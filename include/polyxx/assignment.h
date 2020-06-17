#pragma once

#include <iosfwd>

#include "../assignment.h"
#include "context.h"
#include "value.h"
#include "variable.h"

namespace poly {

  /**
   * Implements a wrapper for lp_assignment_t.
   */
  class Assignment {
    /** The actual assignment. */
    lp_assignment_t mAssignment;

   public:
    /** Construct an empty assignment with a custom context. */
    Assignment(const Context& c);
    /** Construct an empty assignment. */
    Assignment() : Assignment(Context::get_context()) {}
    /** Custom destructor. */
    ~Assignment();

    /** Get a non-const pointer to the internal lp_assignment_t. Handle with
     * care!
     */
    lp_assignment_t* get_internal();
    /** Get a const pointer to the internal lp_assignment_t. */
    const lp_assignment_t* get_internal() const;

    /** Assign variable to the given value. */
    void set(const Variable& var, const Value& value);
    /** Unassign the given variable. */
    void unset(const Variable& var);
    /** Check whether the variable is assigned. */
    bool has(const Variable& var) const;
    /** Retrieve a value from the Assignment. */
    const Value& get(const Variable& var) const;
    /** Clear the assignment. */
    void clear();
  };

  /** Stream the given Assignment to an output stream. */
  std::ostream& operator<<(std::ostream& os, const Assignment& a);

}  // namespace poly
