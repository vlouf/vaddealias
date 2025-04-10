#ifndef CAPPI_H
#define CAPPI_H

#include "pch.h"
using namespace bom;

auto compute_roi(const float range, const float beamwidth, const float max_roi) -> float;
auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t;
auto find_ray(array1<angle> const& rays, angle target) -> size_t;
auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , float roi = 2500.f
    , const float beamwidth = 1.f
    ) -> array2f;

#endif