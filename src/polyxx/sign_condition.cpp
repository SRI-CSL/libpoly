#include "polyxx/sign_condition.h"

#include <cassert>
#include <iostream>

namespace poly {

  std::ostream& operator<<(std::ostream& os, const SignCondition& sc) {
    switch (sc) {
      case SignCondition::LT:
        os << "<";
        break;
      case SignCondition::LE:
        os << "<=";
        break;
      case SignCondition::EQ:
        os << "=";
        break;
      case SignCondition::NE:
        os << "!=";
        break;
      case SignCondition::GT:
        os << ">";
        break;
      case SignCondition::GE:
        os << ">=";
        break;
    }
    return os;
  }

  bool evaluate_sign_condition(SignCondition sc, int sgn) {
    assert(-1 <= sgn && sgn <= 1);
    return lp_sign_condition_consistent(to_sign_condition(sc), sgn);
  }

}  // namespace poly
