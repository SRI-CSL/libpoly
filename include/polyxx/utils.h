#pragma once

#include <functional>
#include <memory>

namespace poly {

  /** Generic type alias for a unique_ptr with a deleter function. */
  template <typename T>
  using deleting_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

}  // namespace poly
