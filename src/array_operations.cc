#include "array_operations.h"

// template <typename T>
// auto argmin(const vector<T>& x) -> size_t{
//   size_t min_index = std::distance(x.begin(), std::min_element(x.begin(), x.end()));
//   return min_index;
// }


auto copy_mask(array2f const& ref, array2f& dest) -> void
{
  for (size_t y = 0; y < ref.extents().y; ++y)
  {
    for (size_t x = 0; x < ref.extents().x; ++x)
    {
      if(std::isnan(ref[y][x]) || std::fabs(ref[y][x] - undetect) < 0.001f)
        dest[y][x] = nodata;
      if(std::fabs(dest[y][x] - undetect) < 0.001f)
        dest[y][x] = nodata;
    }
  }
}

auto flip(array1d& data) -> void
{
  auto copy = array1d{data};
  auto nlen = data.size() - 1;
  for(size_t i = 0; i < copy.size(); i++)
  {
    data[i] = copy[nlen - i];
  }
}

auto flipud(array2f& data) -> void
{
  auto orig = array2f{data};
  auto nlen = data.extents().y - 1;
  for (size_t y = 0; y < data.extents().y; ++y)
  {
    for (size_t x = 0; x < data.extents().x; ++x)
    {
      data[y][x] = orig[nlen - y][x];
    }
  }
}

auto fliplr(array2f& data) -> void
{
  auto orig = array2f{data};
  auto nlen = data.extents().x - 1;
  for (size_t y = 0; y < data.extents().y; ++y)
  {
    for (size_t x = 0; x < data.extents().x; ++x)
    {
      data[y][x] = orig[y][nlen - x];
    }
  }
}

auto mean(const array1f& x) -> float {
  auto sum = 0.0f;
  auto count = 0;

  for (const auto& value : x) {
    if (std::isnan(value))
      continue;

    sum += value;
    count++;
  }

  if (count == 0)
    return nodata;

  return sum / count;
}
