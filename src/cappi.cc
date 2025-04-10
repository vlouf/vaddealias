#include "cappi.h"

auto compute_roi(const float range, const float beamwidth, const float max_roi) -> float{
  float roi = range * beamwidth * (M_PI / 180.0);
  // if(roi > max_roi)
  //   roi = max_roi;
  // else 
  if(roi < 10.f)
    roi = 10.;
  return roi;
}

auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t{
  auto ipos = std::upper_bound(bins.begin(), bins.end(), target, [](auto& l, auto& r) { return l < r.ground_range; });
  if (ipos == bins.end())
    return bins.size();
  if (ipos == bins.begin())
    return 0;
  if (std::fabs(ipos->ground_range - target) < std::fabs((ipos - 1)->ground_range - target))
    return ipos - bins.begin();
  return ipos - bins.begin() - 1;
}

auto find_ray(array1<angle> const& rays, angle target) -> size_t{
  // HACK - not handling 0/360 well
  auto ipos = std::upper_bound(rays.begin(), rays.end(), target);
  if (ipos == rays.end())
    return rays.size() - 1;
  if (ipos == rays.begin())
    return 0;
  if (abs(*ipos - target) < abs(*(ipos - 1) - target))
    return ipos - rays.begin();
  return ipos - rays.begin() - 1;
}

auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , float max_roi
    , const float beamwidth
    ) -> array2f
{
  auto roi = max_roi;
  auto cappi = array2f{latlons.extents()};

  size_t topscan=0;
  size_t botscan=0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++){
    if(vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++){
    if(vol.sweeps[iscan].beam.elevation() < vol.sweeps[botscan].beam.elevation())
      botscan = iscan;
  }

  for (size_t y = 0; y < cappi.extents().y; ++y){
    for (size_t x = 0; x < cappi.extents().x; ++x){
      auto ll = latlons[y][x];
      auto br = wgs84.latlon_to_bearing_range(vol.location, ll);
      br.first = br.first.normalize();

      int lwr_scan = -1, upr_scan = -1;
      int lwr_bin, upr_bin;
      float lwr_val = nodata;
      float upr_val = nodata;
      float lwr_dist = nodata;
      float upr_dist = nodata;
      for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan)
      {
        auto& scan = vol.sweeps[iscan];

        auto ibin = find_ground_range_bin(scan.bins, br.second);
        if (ibin >= scan.bins.size())
          continue;
        roi = compute_roi(scan.bins[ibin].ground_range, beamwidth, max_roi);

        auto alt_dist = scan.bins[ibin].altitude - altitude;
        if (alt_dist > max_alt_diff)
          continue;
        if(iscan == topscan || alt_dist < - max_alt_diff)
          continue;

        auto iray = find_ray(scan.rays, br.first);
        auto val = scan.data[iray][ibin];
        if (std::isnan(val) || std::fabs(val - undetect) < .1f)
          continue;

        if (alt_dist <= 0)
        {
          if (lwr_scan == -1 || scan.bins[ibin].altitude > vol.sweeps[lwr_scan].bins[lwr_bin].altitude)
          {
            lwr_scan = iscan;
            lwr_bin = ibin;
            lwr_val = val;
            lwr_dist = -alt_dist;
          }
        } else {
          if (upr_scan == -1 || scan.bins[ibin].altitude < vol.sweeps[upr_scan].bins[upr_bin].altitude)
          {
            upr_scan = iscan;
            upr_bin = ibin;
            upr_val = val;
            upr_dist = alt_dist;
          }
        }
      }

      if (lwr_scan != -1 && upr_scan != -1)
      {
        if (lwr_dist > 0.0f)
        {
          auto idw_lwr = 1.0 / std::pow(lwr_dist, idw_pwr);
          auto idw_upr = 1.0 / std::pow(upr_dist, idw_pwr);
          auto norm = idw_lwr + idw_upr;
          cappi[y][x] = lwr_val * (idw_lwr / norm) + upr_val * (idw_upr / norm);
        }
        else if (lwr_scan != 0)
          cappi[y][x] = lwr_val;
      }
      else if (lwr_scan != -1 && lwr_dist < roi)
        cappi[y][x] = lwr_val;
      else if (upr_scan != -1 && upr_dist < roi)
        cappi[y][x] = upr_val;
      else
        cappi[y][x] = nodata;
    }
  }

  return cappi;
}
