#include "polyxx/polynomial_utils.h"

namespace poly {

  VariableCollector::VariableCollector() {
    lp_variable_list_construct(&mVarList);
  }
  VariableCollector::~VariableCollector() {
    lp_variable_list_destruct(&mVarList);
  }
  void VariableCollector::operator()(const Polynomial& p) {
    lp_polynomial_get_variables(p.get_internal(), &mVarList);
  }
  std::vector<Variable> VariableCollector::get_variables() const {
    std::vector<Variable> res;
    for (std::size_t i = 0; i < lp_variable_list_size(&mVarList); ++i) {
      res.emplace_back(mVarList.list[i]);
    }
    return res;
  }

}  // namespace poly