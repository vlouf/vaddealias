#ifndef UTIL_IMPORT
#define UTIL_IMPORT

#include "pch.h"
using namespace bom;

// template <typename T>
// auto argmin(const vector<T>& x ) -> size_t;

template <typename T>
auto argmin2(const vector<T>& x, const T val) -> size_t{
  vector<T> x2;
  for(size_t i=0; i<x.size(); i++){
    x2.push_back(std::abs(x[i] - val));
  }
  auto min_it = std::min_element(x2.begin(), x2.end());
  size_t min_index = std::distance(x2.begin(), min_it);
  return min_index;
}

auto copy_mask(array2f const& ref, array2f& dest) -> void;
auto flip(array1d& data) -> void;
auto flipud(array2f& data) -> void;
auto fliplr(array2f& data) -> void;
auto mean(const array1f& x) -> float;

#endif