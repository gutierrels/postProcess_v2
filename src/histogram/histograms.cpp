
#include "histogram/histograms.hh"
#include <cstdio>

HistogramSet::HistogramSet(const SimConfig &cfg, const DetectorGeometry &geo,
                           size_t nPairs)
    : cfg(cfg), geo(geo), nPairs(nPairs) {

  histP.resize(nPairs, 0);
  histR.resize(nPairs, 0);
  histPx.resize(nPairs * cfg.nBinsX * 2, 0);
  histPy.resize(nPairs * cfg.nBinsY * 2, 0);
  histM.resize(geo.totalMods, 0);
  histMx.resize(geo.totalMods * cfg.nBinsX, 0);
  histMy.resize(geo.totalMods * cfg.nBinsY, 0);

  if (cfg.generateHistogram == GENERATE_HISTOGRAM::LOR_INDEX) {
    histLOR.resize(cfg.N * cfg.N * cfg.N * cfg.N * nPairs, 0);
  }

  nBinsE = 300;
  double de = (geo.emax - geo.emin) / static_cast<double>(nBinsE);
  dekeV = de / 1.0e3;
  histMe.resize(geo.totalMods * nBinsE, 0);

  nBinsZ = 10;
  dz = cfg.detectorDepth / static_cast<double>(nBinsZ);
  histMez.resize(nBinsZ * geo.totalMods * nBinsE, 0);

  histNormPairSize =
      cfg.nBinsNormX * cfg.nBinsNormY * cfg.nBinsNormX * cfg.nBinsNormY;
  histNormCentralLOR =
      (cfg.nBinsNormX * cfg.nBinsNormY + 1) *
      (cfg.nBinsNormX * (cfg.nBinsNormY / 2) + cfg.nBinsNormX / 2);
  size_t histNormSize = histNormPairSize * nPairs;
  histNorm.resize(histNormSize, 0.0);
}

void HistogramSet::accumulateCoincidence(
    const single &s1, const single &s2, const std::array<double, 3> &p1,
    const std::array<double, 3> &p2, const coincidence &c, unsigned normBin1X,
    unsigned normBin1Y, unsigned normBin2X,
                             unsigned normBin2Y, size_t lorIdx) {

  ++histP[c.pair];
  ++histPx[2 * c.pair * cfg.nBinsX + c.xPosition1];
  ++histPx[2 * c.pair * cfg.nBinsX + cfg.nBinsX + c.xPosition2];
  ++histPy[2 * c.pair * cfg.nBinsY + c.yPosition1];
  ++histPy[2 * c.pair * cfg.nBinsY + cfg.nBinsY + c.yPosition2];

  const unsigned histPairIndex = c.pair;
  const unsigned pairFirstIndex = histPairIndex * histNormPairSize;
  const unsigned binFirstIndex = (normBin1Y * cfg.nBinsNormX + normBin1X) *
                                 cfg.nBinsNormX * cfg.nBinsNormY;
  histNorm[pairFirstIndex + binFirstIndex + normBin2Y * cfg.nBinsNormX +
           normBin2X] += 1.0;

  if (cfg.generateHistogram == GENERATE_HISTOGRAM::LOR_INDEX) {
    ++histLOR[lorIdx];
  }

  ++histM[s1.module];
  ++histM[s2.module];

  ++histMx[s1.module * cfg.nBinsX + c.xPosition1];
  ++histMx[s2.module * cfg.nBinsX + c.xPosition2];

  ++histMy[s1.module * cfg.nBinsY + c.yPosition1];
  ++histMy[s2.module * cfg.nBinsY + c.yPosition2];

  int ebin = (s1.e - geo.eminkeV) / dekeV;
  if (ebin >= 0 && ebin < static_cast<int>(nBinsE)) {
    ++histMe[s1.module * nBinsE + ebin];

    int zbin = static_cast<int>(p1[2] / dz);
    if (zbin >= 0 && zbin < static_cast<int>(nBinsZ)) {
      ++histMez[s1.module * nBinsE * nBinsZ + ebin * nBinsZ + zbin];
    }
  }

  ebin = (s2.e - geo.eminkeV) / dekeV;
  if (ebin >= 0 && ebin < static_cast<int>(nBinsE)) {
    ++histMe[s2.module * nBinsE + ebin];

    int zbin = static_cast<int>(p2[2] / dz);
    if (zbin >= 0 && zbin < static_cast<int>(nBinsZ)) {
      ++histMez[s2.module * nBinsE * nBinsZ + ebin * nBinsZ + zbin];
    }
  }
}

void HistogramSet::accumulateSingle(const single &s,
                                    const std::array<double, 3> &p) {
  int ebin = (s.e - geo.eminkeV) / dekeV;
  if (ebin >= 0 && ebin < static_cast<int>(nBinsE)) {
    ++histMe[s.module * nBinsE + ebin];

    int zbin = static_cast<int>(p[0] / dz);
    if (zbin >= 0 && zbin < static_cast<int>(nBinsZ)) {
      ++histMez[s.module * nBinsE * nBinsZ + ebin * nBinsZ + zbin];
    }
  }
}
