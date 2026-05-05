
#ifndef POSTPROCESS_HISTOGRAM_HISTOGRAMS_HH
#define POSTPROCESS_HISTOGRAM_HISTOGRAMS_HH

#include "common/types.hh"
#include "config/config.hh"
#include "geometry/geometry.hh"
#include <map>
#include <string>
#include <vector>

struct HistogramSet {
  std::vector<unsigned> histP;
  std::vector<unsigned> histR;
  std::vector<unsigned> histPx;
  std::vector<unsigned> histPy;
  std::vector<unsigned> histM;
  std::vector<unsigned> histMx;
  std::vector<unsigned> histMy;
  std::vector<float> histLOR;

  std::vector<unsigned> histMe;
  std::vector<unsigned> histMez;

  std::vector<float> histNorm;


  const SimConfig &cfg;
  const DetectorGeometry &geo;
  size_t nPairs;

  double dekeV;
  double dz;
  size_t nBinsE;
  size_t nBinsZ;

  size_t histNormPairSize;
  size_t histNormCentralLOR;

  HistogramSet(const SimConfig &cfg, const DetectorGeometry &geo,
               size_t nPairs);

  void accumulateCoincidence(const single &s1, const single &s2,
                             const std::array<double, 3> &p1,
                             const std::array<double, 3> &p2,
                             const coincidence &c, unsigned normBin1X,
                             unsigned normBin1Y, unsigned normBin2X,
                             unsigned normBin2Y, size_t lorIdx);

  void accumulateSingle(const single &s, const std::array<double, 3> &p);

  void writeAll(const std::string &outputDir, double maxTimestamp,
                unsigned long nCoincidences) const;

private:
  void writeModuleHistograms(const std::string &outputDir) const;
};

#endif
