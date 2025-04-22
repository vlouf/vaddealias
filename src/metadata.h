#ifndef METADATA_H
#define METADATA_H

#include "pch.h"
using namespace bom;

struct AttributeValues {
    std::string units;
    std::string standard_name;
    std::string long_name;
};

auto get_azimuth(sweep swp) -> vector<double>;
auto get_date() -> std::string;
auto get_elevation(io::odim::polar_volume const vol_odim) -> array1f;
auto get_lowest_sweep_time(io::odim::polar_volume const vol_odim) -> string;
auto get_nyquist(io::odim::polar_volume const vol_odim) -> array1f;
auto get_range(sweep swp) -> vector<double>;
auto init_altitudes(io::configuration const& config) -> array1f;
auto set_nc_var_attrs(io::nc::variable& varid, const std::string moment) -> void ;
auto str_to_tm(const std::string& datetime) -> std::tm;

#endif