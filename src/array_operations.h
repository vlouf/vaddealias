#ifndef UTIL_IMPORT
#define UTIL_IMPORT

#include "pch.h"
using namespace bom;

auto argmin(array1f const& x) -> int;
auto copy_mask(array2f const& ref, array2f& dest) -> void;
auto flip(array1d& data) -> void;
auto flipud(array2f& data) -> void;
auto fliplr(array2f& data) -> void;
auto mean(const array1f& x) -> float;

#endif