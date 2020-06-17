#include "polyxx/assignment.h"

namespace poly {

  Assignment::Assignment(const Context& c) {
    lp_assignment_construct(get_internal(), c.get_variable_db());
  }
  Assignment::~Assignment() { lp_assignment_destruct(get_internal()); }

  lp_assignment_t* Assignment::get_internal() { return &mAssignment; }
  const lp_assignment_t* Assignment::get_internal() const {
    return &mAssignment;
  }

  void Assignment::set(const Variable& var, const Value& value) {
    lp_assignment_set_value(get_internal(), var.get_internal(),
                            value.get_internal());
  }
  void Assignment::unset(const Variable& var) {
    lp_assignment_set_value(get_internal(), var.get_internal(), nullptr);
  }
  bool Assignment::has(const Variable& var) const {
    return lp_assignment_get_value(get_internal(), var.get_internal())->type !=
           LP_VALUE_NONE;
  }
  const Value& Assignment::get(const Variable& var) const {
    return *detail::cast_from(
        lp_assignment_get_value(get_internal(), var.get_internal()));
  }
  void Assignment::clear() {
    const lp_variable_db_t* var_db = get_internal()->var_db;
    lp_assignment_destruct(get_internal());
    lp_assignment_construct(get_internal(), var_db);
  }

  std::ostream& operator<<(std::ostream& os, const Assignment& a) {
    return stream_ptr(os, lp_assignment_to_string(a.get_internal()));
  }

}  // namespace poly