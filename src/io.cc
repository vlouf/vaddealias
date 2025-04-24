#include "io.h"

auto read_moment(io::odim::polar_volume const vol_odim, string moment, io::configuration const& config) -> volume{
  auto vol = volume{};

  vol.location.lat = vol_odim.latitude() * 1_deg;
  vol.location.lon = vol_odim.longitude() * 1_deg;
  vol.location.alt = vol_odim.height();

  for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan)
  {
    auto scan_odim = vol_odim.scan_open(iscan);
    auto scan = sweep{};

    // if(std::fabs(scan_odim.elevation_angle() - 90) < 0.1f)
    // {
    //   // std::cout << "Skipping 90 deg scan" << std::endl;
    //   continue;
    // }

    scan.beam = radar::beam_propagation{vol.location.alt, scan_odim.elevation_angle() * 1_deg};

    scan.bins.resize(scan_odim.bin_count());
    auto range_scale = scan_odim.range_scale();
    auto range_start = scan_odim.range_start() * 1000 + range_scale * 0.5;
    for (size_t i = 0; i < scan.bins.size(); ++i)
    {
      scan.bins[i].slant_range = range_start + i * range_scale;
      std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) = scan.beam.ground_range_altitude(scan.bins[i].slant_range);
    }

    scan.rays.resize(scan_odim.ray_count());
    auto ray_scale = 360_deg / scan.rays.size();
    auto ray_start = scan_odim.ray_start() * 1_deg + ray_scale * 0.5;
    for (size_t i = 0; i < scan.rays.size(); ++i)
      scan.rays[i] = ray_start + i * ray_scale;

    for (size_t idata = 0; idata < scan_odim.data_count(); ++idata)
    {
      auto data_odim = scan_odim.data_open(idata);
      if (data_odim.quantity() != moment)
        continue;

      scan.data.resize(vec2z{(size_t) scan_odim.bin_count(), (size_t) scan_odim.ray_count()});
      if(moment.compare(config["velocity"]) != 0)
      {
        data_odim.read_unpack(scan.data.data(), nodata, nodata);
      }
      else
      {
        data_odim.read_unpack(scan.data.data(), undetect, nodata);
      }

      vol.sweeps.push_back(std::move(scan));
      break;
    }
  }

  return vol;
}

auto read_global_seamask(string const filename) -> seamask{

  auto dset = io::nc::file{filename, io_mode::read_only};
  size_t nx = dset.lookup_dimension("longitude").size();
  size_t ny = dset.lookup_dimension("latitude").size();

  auto landsea = seamask{vec2z{nx, ny}};
  dset.lookup_variable("latitude").read(landsea.lat);
  dset.lookup_variable("longitude").read(landsea.lon);
  dset.lookup_variable("elevation").read(landsea.mask);

  return landsea;
}

auto check_is_ocean(seamask const& landsea, latlon loc) -> bool{
  vector<float> lon, lat;
  for(auto l: landsea.lon) lon.push_back(l);
  for(auto l: landsea.lat) lat.push_back(l);
  int ilon = argmin2(lon, float(loc.lon.degrees()));
  int ilat = argmin2(lat, float(loc.lat.degrees()));

  return landsea.mask[ilat][ilon] < 0;
}

auto read_refl_corrected(io::odim::polar_volume const vol_odim, io::configuration const& config) -> volume{
  // Read radar file
  auto dbzh = read_moment(vol_odim, "DBZH", config);
  auto dbzh_clean = read_moment(vol_odim, "DBZH_CLEAN", config);

  // Read seamask.
  string filename = config.optional("topography", "/opt/swirl/data/AU_elevation_map.nc");
  auto landsea = read_global_seamask(filename);
  auto radarloc = latlon{dbzh.location.lon, dbzh.location.lat};

  for(size_t i=0; i<dbzh.sweeps.size(); i++){
    auto elev = dbzh.sweeps[i].beam.elevation();
    if(elev > 1_deg * 6)
      continue;

    // std::cout << "Sweep: " << elev << std::endl;
    for(size_t k=0; k<dbzh.sweeps[i].bins.size(); k++){
      if(dbzh.sweeps[i].bins[k].altitude > 4000)
        break;

      for(size_t j=0; j<dbzh.sweeps[i].rays.size(); j++){

        auto r0 = dbzh.sweeps[i].data[j][k];
        if(std::isnan(r0))
          continue;
        if(std::abs(r0 - undetect) < 0.01)
          continue;

        // we may have dbzh, but not dbzh_clean
        if (!dbzh_clean.sweeps.empty()){
          auto r1 = dbzh_clean.sweeps[i].data[j][k];
          if(!std::isnan(r1))
            continue;
          if(std::abs(r1 - undetect) > 0.01)
            continue;
        }

        auto gate_latlon = wgs84.bearing_range_to_latlon(
          radarloc,
          dbzh.sweeps[i].rays[j],
          dbzh.sweeps[i].bins[k].ground_range
        );

        if(check_is_ocean(landsea, gate_latlon))
          dbzh.sweeps[i].data[j][k] = nodata;
      }
    }
  }
  return dbzh;
}


auto read_vad(std::filesystem::path const& filename) -> vadset{
  std::ifstream file(filename);
  std::string line;
  vadset data;

  // Skip the first five lines
  for (int i = 0; i < 5; ++i) {
    std::getline(file, line);
  }

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<float> values((std::istream_iterator<float>(iss)), std::istream_iterator<float>());
    if (values.size() >= 10) {
      data.z.push_back(values[0] * 1000); // Convert to meters
      data.npts.push_back(static_cast<int>(values[1]));
      data.u0.push_back(values[3]);
      data.v0.push_back(values[4]);
      data.w0.push_back(values[5]);
      data.vt.push_back(values[6]);
      data.div.push_back(values[7]);
      data.det.push_back(values[8]);
      data.des.push_back(values[9]);
    }
  }

  // Replace -999.0 with NaN
  std::replace(data.z.begin(), data.z.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.u0.begin(), data.u0.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.v0.begin(), data.v0.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.w0.begin(), data.w0.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.vt.begin(), data.vt.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.div.begin(), data.div.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.det.begin(), data.det.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());
  std::replace(data.des.begin(), data.des.end(), -999.0f, std::numeric_limits<float>::quiet_NaN());

  return data;
}

