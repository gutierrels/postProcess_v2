
#ifndef POSTPROCESS_GEOMETRY_HH
#define POSTPROCESS_GEOMETRY_HH

#include "common/types.hh"
#include "config/config.hh"
#include <array>
#include <vector>

struct DetectorGeometry {
  // Physical dimensions (always needed)
  double physicalDetectorSizeX;
  double physicalDetectorSizeY;
  double physicalRingRad;
  double physicalRingDistance;
  unsigned nPhysicalRings;
  unsigned modPerRing;

  // Working dimensions (logical or physical)
  double detectorSizeX;
  double detectorSizeY;
  double ringRad;
  double ringDistance;
  unsigned nRings;
  unsigned totalMods;
  unsigned nPhysicalTotalMods;

  // Derived bin sizes
  double detectorSizeX05;
  double detectorSizeY05;
  double dBinSizeX;
  double dBinSizeY;
  double dBinSizeNormX;
  double dBinSizeNormY;

  // Logical detector geometry
  double logicalDetectorSizeY;

  // Energy window (eV)
  double emin;
  double emax;
  double eminkeV;

  // Rotation matrices per module and axial offsets
  std::vector<std::array<double, 9>> Rz;
  std::vector<double> Dz;
  std::vector<double> physicalDz;

  static DetectorGeometry build(SimConfig &cfg);
};

#endif
