#ifndef IO_H
#define IO_H

#include "pch.h"
#include "array_operations.h"

using namespace bom;

struct seamask{
    array1f lat;
    array1f lon;
    array2i mask;
    seamask(vec2z shape) : lat{shape.y}, lon{shape.x}, mask{shape} { }
};

auto read_moment(io::odim::polar_volume const vol_odim, string moment, io::configuration const& config) -> volume;
auto read_global_seamask(string const filename) -> seamask;
auto check_is_ocean(seamask const& landsea, latlon loc) -> bool;
auto read_refl_corrected(io::odim::polar_volume const vol_odim, io::configuration const& config) -> volume;
auto read_vad(std::filesystem::path const& filename) -> vadset;

#endif