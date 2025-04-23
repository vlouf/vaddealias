#include <getopt.h>
#include "pch.h"

#include "array_operations.h"
#include "cappi.h"
#include "corrections.h"
#include "metadata.h"
#include "io.h"
#include "brox/brox_optic_flow.h"

using namespace bom;

constexpr auto example_config =
R"(# example layered-flow config

# domain projection
proj4 "+proj=aea +lat_1=-32.2 +lat_2=-35.2 +lon_0=151.209 +lat_0=-33.7008 +a=6378137 +b=6356752.31414 +units=m"

# grid size
size "301 301"

# top left coordinates
left_top "-150500 150500"

# grid resolution
cell_delta "1000 -1000"

# horizontal grid units
units m

# altitude of lowest layer (m)
altitude_base 0.0

# altitude step between layers (m)
altitude_step 500.0

# number of layers
layer_count 41

# radar moment to generate CAPPIs from
moment DBZH

# radar vel moment to generate CAPPIs from
velocity VRADDH

# whether to output the cappis as well as flow fields
output_cappis true

# whether to output the flow magnitude and angle fields
output_polar true

# maximum distance from CAPPI altitude to use reflectivities
max_alt_dist 20000

# exponent for inverse distance weighting when interpolating between vertical levels (2 is a good default)
idw_pwr 2.0

# threshold out cappis to this minimum DBZ before tracking
min_dbz 20

# speckle filter: suppress pixels with less than this many non-zero neighbours (3x3)
speckle_min_neighbours 3

# speckle filter: number of times to apply speckle filter
speckle_iterations 3

# Matrix orientation
origin xy

# parameters for optical flow algorithm
optical_flow
{
  alpha 80
  gamma 7.0
  scales 100
  zfactor 0.5
  tol 0.005
  initer 3
  outiter 15
}

# Land/sea mask file
topography "/opt/swirl/data/AU_elevation_map.nc"

)";

constexpr auto try_again = "try --help for usage instructions\n";
constexpr auto usage_string =
R"(Optical flow tracking of radar volume at multiple altitudes

usage:
  track-layers [options] config.conf lag1.vol.h5 lag0.vol.h5 out.nc

available options:
  -h, --help
      Show this message and exit

  -g, --generate
      Output a sample configuration file and exit

  -t, --trace=level
      Set logging level [log]
        none | status | error | warning | log | debug
)";

constexpr auto short_options = "hgt:";
constexpr struct option long_options[] =
{
    { "help",     no_argument,       0, 'h' }
  , { "generate", no_argument,       0, 'g' }
  , { "trace",    required_argument, 0, 't' }
  , { 0, 0, 0, 0 }
};

auto read_volume(
  const std::filesystem::path& path,
  const io::configuration& config,
  const bool ismainfile
) -> radarset {
  // Logging: Replace with appropriate logging mechanism
  // logDebug("Reading: {}", path.string());
  io::odim::polar_volume vol_odim{path, io_mode::read_only};

  radarset dset;
  const std::string& reflname = config["moment"];
  std::string topo_fname = config.optional("topography", "");

  // Get Nyquist
  dset.elevation = get_elevation(vol_odim);
  dset.nyquist = get_nyquist(vol_odim);
  dset.lowest_sweep_time = get_lowest_sweep_time(vol_odim);

  const auto& velocity_moment = config["velocity"];
  dset.vradh = read_moment(vol_odim, velocity_moment, config);

  if (topo_fname.empty()) {
    std::cout << "No topgraphy provided. Not correcting for sea-clutter" << std::endl;
    dset.dbzh = read_moment(vol_odim, reflname, config);
  } else {
    if (ismainfile && reflname == "DBZH") {
      // Correct for sea-clutter
      dset.dbzh = read_refl_corrected(vol_odim, config);
    } else {
      // Read reflectivity without any modification.
      dset.dbzh = read_moment(vol_odim, reflname, config);
    }
  }

  const auto& attributes = vol_odim.attributes();
  dset.source = attributes["source"].get_string();
  dset.date = attributes["date"].get_string();
  dset.time = attributes["time"].get_string();
  dset.beamwidth = attributes["beamwH"].get_real();

  // Logging: Replace with appropriate logging mechanism
  // logDebug("{} read", path.string());
  return dset;
}

template <typename T>
auto meshgrid(const std::vector<T>& x, const std::vector<T>& y) ->
std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>>{
  // Not used yet.

  std::vector<std::vector<T>> X(y.size(), std::vector<T>(x.size()));
  std::vector<std::vector<T>> Y(y.size(), std::vector<T>(x.size()));

  for (size_t i = 0; i < y.size(); ++i) {
    for (size_t j = 0; j < x.size(); ++j) {
      X[i][j] = x[j];
      Y[i][j] = y[i];
    }
  }
  return {X, Y};
}

template <typename T>
auto argmin(const vector<T>& x, T val) -> size_t{
  vector<T> x2;
  for(size_t i=0; i<x.size(); i++){
    x2.push_back(std::abs(x[i] - val));
  }
  auto min_it = std::min_element(x2.begin(), x2.end());
  size_t min_index = std::distance(x2.begin(), min_it);
  return min_index;
}

auto generate_vad_field(radarset dset2, vadset df) -> vector<array2f>{
  vector<array2f>  vadfield;

  // Create the VAD velocity field.
  for(size_t k=0; k < dset2.vradh.sweeps.size(); k++){
    auto r = get_range(dset2.vradh.sweeps[k]);
    auto azi = get_azimuth(dset2.vradh.sweeps[k]);
    auto el = dset2.vradh.sweeps[k].beam.elevation();
    auto cel = cos(M_PI / 180. * el);
    auto sel = sin(M_PI / 180. * el);

    // Generated anew for every sweep
    auto vrz = array2f{vec2z{r.size(), azi.size()}};  // vec2 dimensions are reversed.
    for(size_t i=0; i<r.size(); i++){
      auto alti = dset2.vradh.sweeps[k].bins[i].altitude;
      auto pos = argmin(df.z, alti);

      for(size_t j=0; j<azi.size(); j++){
        vrz[j][i] = (
          0.5 * r[i] * cel * df.div[pos]
          - df.vt[pos] * sel
          + df.u0[pos] * sin(M_PI / 180. * azi[j]) * cel
          + df.v0[pos] * cos(M_PI / 180. * azi[j]) * cel
          - 0.5 * r[i] * cel * cos(2 * M_PI / 180. * azi[j]) * df.det[pos]
          + 0.5 * r[i] * cel * sin(2 * M_PI / 180. * azi[j]) * df.des[pos]
        );
      }
    }
    vadfield.push_back(vrz);
  }
  return vadfield;
}

auto process_file(
  io::configuration const& config,
  std::filesystem::path const& vad_file,
  std::filesystem::path const& odim_file1,
  std::filesystem::path const& odim_file2
) -> void{
  auto dset1 = read_volume(odim_file1, config, true);
  auto dset2 = read_volume(odim_file2, config, false);
  auto df = read_vad(vad_file);
  auto vadfield = generate_vad_field(dset2, df);
}

auto check_configuration_file(io::configuration const& config) -> bool
{
  bool result = false;

  if(config["origin"].string().compare("ij") == 0 || config["origin"].string().compare("xy") == 0)
  {
    result = true;
  }
  else
  {
    std::cout << "Invalid parameter in configuration file. Origin = 'ij' | 'xy'. You gave: " << config["origin"].string() << std::endl;
    throw std::invalid_argument("Invalid values for parameter in configuration file ('origin')");
  }

  return result;
}

int main(int argc, char* argv[])
{
  try
  {
    // process command line
    while (true)
    {
      int option_index = 0;
      int c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
        break;
      switch (c)
      {
      case 'h':
        std::cout << usage_string;
        return EXIT_SUCCESS;
      case 'g':
        std::cout << example_config;
        return EXIT_SUCCESS;
      case 't':
        trace::set_min_level(from_string<trace::level>(optarg));
        break;
      case '?':
        std::cerr << try_again;
        return EXIT_FAILURE;
      }
    }

    if (argc - optind != 4)
    {
      std::cerr << "missing required parameter\n" << try_again;
      return EXIT_FAILURE;
    }

    auto config = io::configuration{std::ifstream{argv[optind+0]}};
    if(check_configuration_file(config) != true)
      return EXIT_FAILURE;

    process_file(
          config
        , argv[optind+1]
        , argv[optind+2]
        , argv[optind+3]
        );
  }
  catch (std::exception& err)
  {
    trace::error("fatal exception: {}", format_exception(err));
    return EXIT_FAILURE;
  }
  catch (...)
  {
    trace::error("fatal exception: (unknown exception)");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
