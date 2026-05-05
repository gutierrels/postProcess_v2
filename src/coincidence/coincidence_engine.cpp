
#include "coincidence/coincidence_engine.hh"
#include <cmath>

CoincidenceEngine::CoincidenceEngine(unsigned method,
                                     const std::vector<detPair> &pairs,
                                     unsigned modPerRing, bool logicalMode)
    : coinMethod(method), pairIndexes(pairs), modPerRing(modPerRing),
      useLogicalDetectors(logicalMode) {}

CoincidenceResult
CoincidenceEngine::findCoincidence(const std::vector<single> &window,
                                   size_t windowSize) const {
  CoincidenceResult result = {0, static_cast<unsigned>(pairIndexes.size())};

  if (windowSize <= 1) {
    return result;
  }

  // Get first single history number
  unsigned long long hist = window[0].hist;

  // Find the closest single scored by a related detector
  double emaxSingle = 0.0;
  double d511 = 1.0e35;

  unsigned nextSingDevMod =
      useLogicalDetectors ? window[0].module
                          : sim2devModule(modPerRing, window[0].module);
  unsigned localPair;

  switch (coinMethod) {
  case COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(
          pairIndexes, nextSingDevMod,
          useLogicalDetectors ? window[i].module
                              : sim2devModule(modPerRing, window[i].module));
      if (localPair >= pairIndexes.size())
        continue;

      // Check if can be considered as the next coincidence
      if (emaxSingle < window[i].e) {
        emaxSingle = window[i].e;
        result.iCoincidence = i;
        result.iPair = localPair;
      }
    }
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(
          pairIndexes, nextSingDevMod,
          useLogicalDetectors ? window[i].module
                              : sim2devModule(modPerRing, window[i].module));
      if (localPair >= pairIndexes.size())
        continue;

      // Check if can be considered as the next coincidence
      if (hist == window[i].hist) {
        result.iCoincidence = i;
        result.iPair = localPair;
        break;
      }
    }
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(
          pairIndexes, nextSingDevMod,
          useLogicalDetectors ? window[i].module
                              : sim2devModule(modPerRing, window[i].module));
      if (localPair >= pairIndexes.size())
        continue;

      // Check if can be considered as the next coincidence
      double locald511 = std::fabs(window[i].e - 511.0);
      if (locald511 < d511 && hist == window[i].hist) {
        d511 = locald511;
        result.iCoincidence = i;
        result.iPair = localPair;
      }
    }
    break;
  case COINCIDENCE_METHOD::TAKE_CLOSEST:
  default:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(
          pairIndexes, nextSingDevMod,
          useLogicalDetectors ? window[i].module
                              : sim2devModule(modPerRing, window[i].module));
      if (localPair >= pairIndexes.size())
        continue;

      result.iCoincidence = i;
      result.iPair = localPair;
      break;
    }
    break;
  }

  return result;
}
