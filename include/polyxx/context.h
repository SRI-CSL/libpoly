#pragma once

#include "../polynomial_context.h"
#include "../variable_db.h"
#include "../variable_order.h"
#include "utils.h"

namespace poly {
  /**
   * The wrapper class for the polynomial context.
   * Only supports the lp_Z ring right now.
   */
  class Context {
    /** A variable database. */
    deleting_unique_ptr<lp_variable_db_t> mVariableDB;
    /** A variable ordering. */
    deleting_unique_ptr<lp_variable_order_t> mVariableOrder;
    /** A polynomial context. */
    deleting_unique_ptr<lp_polynomial_context_t> mPolynomialContext;

   public:
    /** Constructs a new (empty) context */
    Context();

    /** Wraps the lp context. */
    explicit Context(lp_polynomial_context_t* ctx);

    /** Copy constructor. */
    Context(const Context& other);

    /** No assignment as this would mess with internal reference counting. */
    Context& operator=(const Context& other) = delete;

    /** Get a non-const pointer to the internal lp_variable_db_t.
     * Handle with care!
     */
    lp_variable_db_t* get_variable_db() const;

    /** Get a non-const pointer to the internal lp_variable_order_t.
     * Handle with care!
     */
    lp_variable_order_t* get_variable_order() const;

    /** Get a non-const pointer to the internal lp_polynomial_context_t.
     * Handle with care!
     */
    lp_polynomial_context_t* get_polynomial_context() const;

    static const Context& get_context();
  };

}  // namespace poly