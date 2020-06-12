#pragma once

#include <iostream>

#include "../variable_db.h"
#include "context.h"

namespace poly {

  /**
   * Implements a wrapper for lp_variable_t from libpoly.
   */
  class Variable {
    /** The actual variable. */
    lp_variable_t mVariable;

   public:
    /** Construct with a null variable. */
    Variable();
    /** Construct from a lp_variable_t. */
    Variable(lp_variable_t var);
    /** Construct a new variable with the given name in the specified context.
     */
    Variable(const Context& c, const char* name);
    /** Construct a new variable with the given name in the default context. */
    Variable(const char* name);

    /** Get the internal lp_variable_t. Note that it's only a type alias for
     * long.
     */
    lp_variable_t get_internal() const;
  };

  namespace detail {
    /** Utility struct for easy printing of variables in custom contexts. */
    struct variable_printer {
      const lp_variable_db_t* var_db;
      lp_variable_t var;
    };
    /** Actually stream a variable_printer object. */
    inline std::ostream& operator<<(std::ostream& os,
                                    const variable_printer& v) {
      return os << lp_variable_db_get_name(v.var_db, v.var);
    }
  }  // namespace detail

  /** Stream the given Variable to an output stream from the default context. */
  std::ostream& operator<<(std::ostream& os, const Variable& v);

  /** Stream the given Variable from the given context.
   * Use as follows: os << stream_variable(c, v) << ...
   */
  inline detail::variable_printer stream_variable(const Context& c,
                                                  const Variable& v);

  /** Compare two variables for equality. */
  bool operator==(const Variable& lhs, const Variable& rhs);
  /** Compare two variables for inequality. */
  bool operator!=(const Variable& lhs, const Variable& rhs);

}  // namespace poly
