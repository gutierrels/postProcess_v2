
#include "geometry/geometry.hh"
#include <cmath>
#include <cstdio>

DetectorGeometry DetectorGeometry::build(SimConfig &cfg) {

  DetectorGeometry geo;

  // Physical dimensions (always needed)
  geo.physicalDetectorSizeX = cfg.header.detectorSizeX;
  geo.physicalDetectorSizeY = cfg.header.detectorSizeY;
  geo.physicalRingRad = cfg.header.detectorDistance / 2.0;
  geo.physicalRingDistance = cfg.header.ringDistance;
  geo.nPhysicalRings = cfg.header.ringNumber;
  geo.modPerRing = cfg.header.moduleNumber;

  // Working dimensions default to physical
  geo.detectorSizeX = geo.physicalDetectorSizeX;
  geo.detectorSizeY = geo.physicalDetectorSizeY;
  geo.ringRad = geo.physicalRingRad;
  geo.ringDistance = geo.physicalRingDistance;
  geo.nRings = geo.nPhysicalRings;

  // Logical detector geometry
  geo.logicalDetectorSizeY = 0.0;
  if (cfg.useLogicalDetectors) {
    const double scannerLength =
        geo.physicalDetectorSizeY * geo.nPhysicalRings +
        geo.physicalRingDistance * geo.nPhysicalRings;

    printf("Scanner length: %f cm (Contando medio gap arriba y medio abajo)\n",
           scannerLength);

    geo.logicalDetectorSizeY =
        scannerLength / static_cast<double>(cfg.nLogicalRings);
    const double logicalDetectorSizeX = geo.logicalDetectorSizeY;

    printf("Physical detector size X: %f cm\n", geo.physicalDetectorSizeX);
    printf("Physical detector size Y: %f cm\n", geo.physicalDetectorSizeY);
    printf("Physical detector rings: %u\n", geo.nPhysicalRings);
    printf("-------------------------------------\n");
    printf("Logical detector size X: %f cm\n", logicalDetectorSizeX);
    printf("Logical detector size Y: %f cm\n", geo.logicalDetectorSizeY);
    printf("Logical detector rings: %u\n", cfg.nLogicalRings);
    printf("-------------------------------------\n");

    // Override working dimensions
    geo.detectorSizeX = logicalDetectorSizeX;
    geo.detectorSizeY = geo.logicalDetectorSizeY;
    geo.ringRad = geo.physicalRingRad;
    geo.ringDistance = 0.0;
    geo.nRings = cfg.nLogicalRings;

    // Override header for LM output
    cfg.header.detectorSizeX = logicalDetectorSizeX;
    cfg.header.detectorSizeY = geo.logicalDetectorSizeY;
    cfg.header.ringDistance = 0.0;
    cfg.header.ringNumber = cfg.nLogicalRings;
  }

  geo.detectorSizeX05 = geo.detectorSizeX / 2.0;
  geo.detectorSizeY05 = geo.detectorSizeY / 2.0;
  geo.dBinSizeX = geo.detectorSizeX / static_cast<double>(cfg.nBinsX);
  geo.dBinSizeY = geo.detectorSizeY / static_cast<double>(cfg.nBinsY);
  geo.dBinSizeNormX =
      geo.detectorSizeX / static_cast<double>(cfg.nBinsNormX);
  geo.dBinSizeNormY =
      geo.detectorSizeY / static_cast<double>(cfg.nBinsNormY);

  // Convert header distances to mm (for LM output)
  cfg.header.detectorSizeX *= 10.0;
  cfg.header.detectorSizeY *= 10.0;
  cfg.header.detectorDistance *= 10.0;
  cfg.header.ringDistance *= 10.0;

  geo.totalMods = geo.modPerRing * geo.nRings;
  geo.nPhysicalTotalMods = geo.modPerRing * geo.nPhysicalRings;

  // Energy window (eV)
  geo.emin = 511.0e3 * (1.0 - cfg.energyWindow);
  geo.emax = 511.0e3 * (1.0 + cfg.energyWindow);
  geo.eminkeV = geo.emin / 1.0e3;

  // Physical Dz (always needed, used for toLogical in logical mode)
  geo.physicalDz.resize(geo.nPhysicalRings);
  for (unsigned iRing = 0; iRing < geo.nPhysicalRings; ++iRing) {
    geo.physicalDz[iRing] =
        iRing * (geo.physicalRingDistance + geo.physicalDetectorSizeY);
  }

  // Working Dz (logical or physical depending on mode)
  geo.Dz.resize(geo.nRings);
  if (cfg.useLogicalDetectors) {
    for (unsigned iRing = 0; iRing < geo.nRings; ++iRing) {
      geo.Dz[iRing] = iRing * (geo.ringDistance + geo.detectorSizeY);
    }
  } else {
    geo.Dz = geo.physicalDz;
  }

  // Rotation matrices per module
  geo.Rz.resize(geo.modPerRing);
  constexpr double pi = 3.14159265359;
  const double angleStep =
      2.0 * pi / static_cast<double>(geo.modPerRing);
  for (unsigned im = 0; im < geo.modPerRing; ++im) {
    double angle = static_cast<double>(im) * angleStep;
    double cangle = cos(angle);
    double sangle = sin(angle);
    geo.Rz[im] = {cangle, -sangle, 0.0, sangle, cangle, 0.0, 0.0, 0.0, 1.0};
  }

  return geo;
}
