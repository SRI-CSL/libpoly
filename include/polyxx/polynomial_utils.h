#pragma once

#include "polynomial.h"
#include "../variable_list.h"

#include <vector>

namespace poly {

  /** Utility class to collect all variables from a sequence of polynomials.
   */
  class VariableCollector {
   private:
    /** Internal variable list. */
    lp_variable_list_t mVarList;

   public:
    /** Default constructor, create the variable list. */
    VariableCollector();
    /** Destructs the variable list. */
    ~VariableCollector();
    /** Adds the variables of the given polynomial to the list. */
    void operator()(const Polynomial& p);
    /** Returns the list of variables as a vector of Variables. */
    std::vector<Variable> get_variables() const;
  };

}  // namespace poly
