
#ifndef POSTPROCESS_COINCIDENCE_ENGINE_HH
#define POSTPROCESS_COINCIDENCE_ENGINE_HH

#include "common/types.hh"
#include <vector>

struct CoincidenceResult {
  unsigned iCoincidence; // 0 = no CoincidenceEvent
  unsigned iPair;
  bool isTriple;
};

class CoincidenceEngine {
public:
  CoincidenceEngine(CoincidenceMethod method,
                    const std::vector<DetectorPair> *pairs, unsigned modPerRing,
                    bool useLogicalDetectors, double emin511keV = 0.0,
                    double emax511KeV = 0.0, double promptGammaMin = 0.0,
                    double promptGammaMax = 0.0);

  CoincidenceResult findCoincidence(const std::vector<SingleEvent> &window,
                                    size_t windowSize) const;

private:
  CoincidenceMethod coinMethod;
  const std::vector<DetectorPair> *pairIndexes;
  unsigned modPerRing;
  bool useLogicalDetectors;
  double emin511KeV;
  double emax511KeV;
  double promptGammaMin;
  double promptGammaMax;

  inline unsigned resolveModule(unsigned rawModule) const {
    return useLogicalDetectors ? rawModule
                               : sim2devModule(modPerRing, rawModule);
  }
};

#endif
