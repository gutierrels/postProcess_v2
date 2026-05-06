
#ifndef POSTPROCESS_CONFIG_HH
#define POSTPROCESS_CONFIG_HH

#include "common/types.hh"
#include <string>

struct SimConfig {
  LMHeader header;

  double eRes;              // fractional (already divided by 100)
  double pRes;              // cm
  double tRes;              // s
  double coincidenceWindow; // s
  double energyWindow;      // fractional (already divided by 100)

  unsigned nBinsX;
  unsigned nBinsY;
  unsigned nBinsNormX;
  unsigned nBinsNormY;
  unsigned N; // = nBinsNormX = nBinsNormY (must be equal)

  double detectorDepth; // cm

  unsigned long long nCounts;
  unsigned long long maxCounts;

  unsigned outputFormat;
  unsigned projectionMethod;
  unsigned coinMethod;

  unsigned generateHistogram;


  std::string pairListFilename;

  bool useLogicalDetectors;
  unsigned nLogicalRings;
  unsigned discardMultiples;
};

SimConfig parseConfig(const std::string &filename);

#endif
