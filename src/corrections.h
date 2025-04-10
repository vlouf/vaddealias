#ifndef CORR_H
#define CORR_H

#include "pch.h"
#include "array_operations.h"
using namespace bom;

auto correct_undetect(volume& vol) -> void;
auto unfold(float v, float vref, float vnq, float vshift) -> float;
auto mad_filter(
      volume &velocity
      , const array1f& nyquist
      , float delta_vmax = 5.f
      , size_t nfilter = 10
    ) -> void;
auto speckle_filter(array2f& data, float min_dbz, int min_neighbours) -> void;

#endif