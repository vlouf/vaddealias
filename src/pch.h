#ifndef PCH_H
#define PCH_H

#include <bom/io/configuration.h>
#include <bom/io/cf.h>
#include <bom/io/nc.h>
#include <bom/io/odim.h>
#include <bom/radar/beam_propagation.h>
#include <bom/array2.h>
#include <bom/ellipsoid.h>
#include <bom/grid_coordinates.h>
#include <bom/grid_transform.h>
#include <bom/map_projection.h>
#include <bom/trace.h>

#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace bom;

constexpr float nodata = std::numeric_limits<float>::quiet_NaN();
constexpr float undetect = -32.0f;

struct bin_info{
  float slant_range;
  float ground_range;
  float altitude;
};

struct sweep{
  radar::beam_propagation beam;
  array1<bin_info>        bins; // @ bin centers
  array1<angle>           rays; // @ ray centers
  array2f                 data;  
};

struct volume{
  latlonalt     location;
  vector<sweep> sweeps;
};

struct radarset{
  volume vradh;
  volume dbzh;
  array1f nyquist;
  array1f elevation;
  string source;
  string date;
  string time;
  string lowest_sweep_time;
  float beamwidth;
};

struct vadset{
  vector<float> z;
  vector<int> npts;
  vector<float> u0;
  vector<float> v0;
  vector<float> w0;
  vector<float> vt;
  vector<float> div;
  vector<float> det;
  vector<float> des;
};

#endif // PCH_H