#include "metadata.h"

auto get_azimuth(sweep swp) -> vector<double>{
  vector<double> azi;
  for(const auto& bin : swp.rays) {
      azi.push_back(bin.degrees());
  }
  return azi;
}

auto get_date() -> std::string {
  std::time_t rawtime = std::time(nullptr);
  std::tm* timeinfo = std::localtime(&rawtime);

  std::stringstream ss;
  ss << std::put_time(timeinfo, "%FT%X");

  return ss.str();
}

auto get_elevation(io::odim::polar_volume const vol_odim) -> array1f{
  const size_t nelev = vol_odim.scan_count();
  auto elevation = array1f{nelev};

  for (size_t iscan = 0; iscan < nelev; ++iscan){
        auto scan_odim = vol_odim.scan_open(iscan);
        elevation[iscan] = scan_odim.elevation_angle();
    }

  return elevation;
}

auto get_lowest_sweep_time(io::odim::polar_volume const vol_odim) -> bom::string{
  const size_t nelev = vol_odim.scan_count();
  auto elev = 90.f;
  bom::timestamp stdate, eddate;

  for (size_t iscan = nelev - 1; iscan > 0; --iscan){
      auto scan_odim = vol_odim.scan_open(iscan);
      if(scan_odim.elevation_angle() < elev){
        stdate = bom::io::odim::strings_to_time(
            scan_odim.attributes()["startdate"].get_string()
          , scan_odim.attributes()["starttime"].get_string()
          );
        eddate = bom::io::odim::strings_to_time(
            scan_odim.attributes()["enddate"].get_string()
          , scan_odim.attributes()["endtime"].get_string()
          );
      }
      elev = scan_odim.elevation_angle();
  }
  auto midpoint_time = stdate + (eddate - stdate) / 2;
  return to_string(midpoint_time);
}

auto get_nyquist(io::odim::polar_volume const vol_odim) -> array1f{
  const size_t nelev = vol_odim.scan_count();
  auto nyquist = array1f{nelev};

  for (size_t iscan = 0; iscan < nelev; ++iscan){
    auto scan_odim = vol_odim.scan_open(iscan);
    if (auto iatt = scan_odim.attributes().find("NI"); iatt != scan_odim.attributes().end())
      nyquist[iscan] = scan_odim.attributes()["NI"].get_real();
    else
      nyquist[iscan] = -9999.;
    }

  return nyquist;
}

auto get_range(sweep swp) -> vector<double>{
  vector<double> r;
  for(const auto& bin : swp.bins) {
      r.push_back(double(bin.slant_range));
  }
  return r;
}

auto init_altitudes(io::configuration const& config) -> array1f
{
  auto alts = array1f{config["layer_count"]};
  auto base = float(config["altitude_base"]);
  auto step = float(config["altitude_step"]);
  for (size_t i = 0; i < alts.size(); ++i)
    alts[i] = base + step * i;
  return alts;
}

auto set_nc_var_attrs(io::nc::variable& varid, const std::string moment) -> void {
  // key, {units, standard_name, long_name}
  static const std::unordered_map<std::string, AttributeValues> attributes = {
    {"latitude", {"degrees_north", "latitude", "latitude_degrees_north"}},
    {"longitude", {"degrees_east", "longitude", "longitude_degrees_east"}},
    {"altitude", {"m", "projection_z_coordinate", ""}},
    {"nyquist", {"m s-1", "", ""}},
    {"reflectivity", {"dBZ", "equivalent_reflectivity_factor", ""}},
    {"velocity", {"m s-1", "", "Doppler velocity"}},
    {"radar_latitude", {"degrees_north", "", "Origin latitude of the radar."}},
    {"radar_longitude", {"degrees_east", "", "Origin longitude of the radar."}},
    {"radar_altitude", {"m", "", "Origin altitude of the radar."}}
  };

  const auto it = attributes.find(moment);
  if (it != attributes.end()) {
    const AttributeValues& attrValues = it->second;
    if (!attrValues.units.empty())
      varid.att_set("units", attrValues.units);
    if (!attrValues.standard_name.empty())
      varid.att_set("standard_name", attrValues.standard_name);
    if (!attrValues.long_name.empty())
      varid.att_set("long_name", attrValues.long_name);
  }
}
