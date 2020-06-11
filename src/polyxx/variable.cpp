#include "polyxx/variable.h"

namespace poly {

  Variable::Variable() : mVariable(lp_variable_null) {}
  Variable::Variable(lp_variable_t var) : mVariable(var) {}
  Variable::Variable(const Context& c, const char* name)
      : mVariable(lp_variable_db_new_variable(c.get_variable_db(), name)) {}
  Variable::Variable(const char* name)
      : Variable(Context::get_context(), name) {}

  lp_variable_t Variable::get_internal() const { return mVariable; }

  std::ostream& operator<<(std::ostream& os, const Variable& v) {
    return os << lp_variable_db_get_name(
               Context::get_context().get_variable_db(), v.get_internal());
  }

  detail::variable_printer stream_variable(const Context& c,
                                           const Variable& v) {
    return detail::variable_printer{c.get_variable_db(), v.get_internal()};
  }

  bool operator==(const Variable& lhs, const Variable& rhs) {
    return lhs.get_internal() == rhs.get_internal();
  }
  bool operator!=(const Variable& lhs, const Variable& rhs) {
    return lhs.get_internal() != rhs.get_internal();
  }

}  // namespace poly
