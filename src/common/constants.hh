
#ifndef POSTPROCESS_COMMON_CONSTANTS_HH
#define POSTPROCESS_COMMON_CONSTANTS_HH

#include <cmath>

namespace constants {

/// Sentinel value for uninitialized timestamps (effectively +∞)
constexpr double SENTINEL_TIME = 1.0e35;

/// FWHM to sigma conversion factor: 2 * sqrt(2 * ln(2))
constexpr double FWHM_TO_SIGMA = 2.355;

/// Conversion factor: 1 keV = 1000 eV
constexpr double EV_PER_KEV = 1.0e3;

/// Conversion factor: 1 Bq ≈ 2.7027e-8 µCi
constexpr double BQ_TO_MICROCURIE = 2.7027027027027e-8;

/// π
constexpr double PI = 3.14159265358979323846;

} // namespace constants

#endif
