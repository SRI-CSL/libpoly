#include "polyxx/utils.h"

#include <cstdlib>
#include <iostream>

namespace poly {
  std::ostream& stream_ptr(std::ostream& os, char* ptr) {
    os << ptr;
    free(ptr);
    return os;
  }
}