#pragma once

#include <functional>
#include <iosfwd>
#include <memory>

namespace poly {

  /** Generic type alias for a unique_ptr with a deleter function. */
  template <typename T>
  using deleting_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

  /** Writes a char pointer to an output stream and frees it afterwards. */
  std::ostream& stream_ptr(std::ostream& os, char* ptr);

}  // namespace poly
