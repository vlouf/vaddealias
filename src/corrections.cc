#include "corrections.h"

auto correct_undetect(volume& vol) -> void {
  for (auto& scan : vol.sweeps) {
    for (size_t ibin = 0; ibin < scan.bins.size(); ++ibin) {
      for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
        if (std::fabs(scan.data[iray][ibin] - undetect) < 0.001f) {
          scan.data[iray][ibin] = nodata;
        }
      }
    }
  }
}

auto unfold(float v, float vref, float vnq, float vshift) -> float
{
  auto delv = v - vref;
  float unfld = 0.0;
  if(std::fabs(delv) < vnq){
    unfld = v;
  } else {
    auto sign = (delv > 0 ? 1 : -1);
    unfld = v - static_cast<int>(delv + sign * vnq / vshift) * vshift;
  }

  return unfld;
}

auto mad_filter(
      volume &velocity
      , const array1f& nyquist
      , float delta_vmax
      , size_t nfilter
    ) -> void
{
  // Mean Average Deviation filtering technique to remove noise from Doppler.
  int cnt = 0;
  for (size_t iscan = 0; iscan < velocity.sweeps.size(); ++iscan)
  {
    auto vnyquist = nyquist[iscan];
    auto vshift = 2 * vnyquist;
    auto& scan = velocity.sweeps[iscan];
    for(size_t iray = 0; iray < scan.rays.size(); ++iray)
    {
      for(size_t ibin = 0; ibin < scan.bins.size() - nfilter; ++ibin)
      {
        auto n1 = ibin;
        auto n2 = (n1 + nfilter < scan.bins.size() ? n1 + nfilter : scan.bins.size());
        array1f vselected{n2 - n1}, vplus{n2 - n1}, vminus{n2 - n1};
        vselected.fill(nodata);
        vplus.fill(nodata);
        vminus.fill(nodata);

        for(size_t i = 0 ; i < n2 - n1 ; ++i)
        {
          vselected[i] = scan.data[iray][i + n1];
          if(vselected[i] >= 0)
            vplus[i] = scan.data[iray][i + n1];
          else
            vminus[i] = scan.data[iray][i + n1];
        }

        auto vmean = mean(vselected);
        if(std::isnan(vmean))
          continue;
        auto vmean_plus = mean(vplus);
        auto vmean_minus = mean(vminus);

        for(size_t i = 0 ; i < n2 - n1 ; ++i)
        {
          if(std::isnan(vselected[i]))
            continue;
          auto vk = vselected[i];
          if(std::fabs(vk - vmean) >= delta_vmax)
          {
            auto vk_unfold = nodata;
            auto dvk = 0.f;
            if(vmean >= 0)
            {
              vk_unfold = unfold(vk, vmean_plus, vnyquist, vshift);
              dvk = std::fabs(vk - vmean_plus);
            }
            else
            {
              vk_unfold = unfold(vk, vmean_minus, vnyquist, vshift);
              dvk = std::fabs(vk - vmean_minus);
            }
            if(std::isnan(vk_unfold))
              continue;
            if(std::fabs(vk_unfold - vmean) < delta_vmax || dvk < delta_vmax)
            {
              cnt += 1;
              velocity.sweeps[iscan].data[iray][ibin + i] = vk_unfold;
            }
            else
            {
              velocity.sweeps[iscan].data[iray][ibin + i] = nodata;
            }
          }
        }
      }
    }
  }
}


auto speckle_filter(array2f& data, float min_dbz, int min_neighbours) -> void
{
  auto copy = data;
  for (size_t y = 1; y < data.extents().y - 1; ++y)
  {
    for (size_t x = 1; x < data.extents().x - 1; ++x)
    {
      if (copy[y][x] <= min_dbz)
        continue;
      auto count
        = (copy[y-1][x-1] > min_dbz ? 1 : 0)
        + (copy[y-1][ x ] > min_dbz ? 1 : 0)
        + (copy[y-1][x+1] > min_dbz ? 1 : 0)
        + (copy[ y ][x-1] > min_dbz ? 1 : 0)
        + (copy[ y ][x+1] > min_dbz ? 1 : 0)
        + (copy[y+1][x-1] > min_dbz ? 1 : 0)
        + (copy[y+1][ x ] > min_dbz ? 1 : 0)
        + (copy[y+1][x+1] > min_dbz ? 1 : 0)
        ;
      if (count < min_neighbours)
        data[y][x] = min_dbz;
    }
  }
}