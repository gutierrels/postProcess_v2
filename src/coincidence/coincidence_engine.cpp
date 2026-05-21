
#include "coincidence/coincidence_engine.hh"
#include "common/constants.hh"
#include <cmath>

CoincidenceEngine::CoincidenceEngine(CoincidenceMethod method,
                                     const std::vector<DetectorPair> *pairs,
                                     unsigned modPerRing, bool logicalMode)
    : coinMethod(method), pairIndexes(pairs), modPerRing(modPerRing),
      useLogicalDetectors(logicalMode) {}

CoincidenceResult
CoincidenceEngine::findCoincidence(const std::vector<SingleEvent> &window,
                                   size_t windowSize) const {
  CoincidenceResult result = {0, static_cast<unsigned>(pairIndexes->size())};

  if (windowSize <= 1) {
    return result;
  }

  // Get first SingleEvent history number
  unsigned long long hist = window[0].hist;

  // Find the closest SingleEvent scored by a related detector
  double emaxSingle = 0.0;
  double d511 = constants::SENTINEL_TIME;

  unsigned nextSingDevMod = resolveModule(window[0].module);
  unsigned localPair;

  switch (coinMethod) {
  case CoincidenceMethod::TAKE_WINNER_IF_ALL_ARE_GOODS: {
    bool allGood = true;
    for (size_t i = 1; i < windowSize; ++i) {
      localPair = getPair(*pairIndexes, nextSingDevMod,
                          resolveModule(window[i].module));
      if (localPair >= pairIndexes->size()) {
        allGood = false;
        break;
      }
    }

    if (allGood) {
      for (size_t i = 1; i < windowSize; ++i) {
        localPair = getPair(*pairIndexes, nextSingDevMod,
                            resolveModule(window[i].module));
        if (emaxSingle < window[i].e) {
          emaxSingle = window[i].e;
          result.iCoincidence = i;
          result.iPair = localPair;
        }
      }
    }
    break;
  }
  case CoincidenceMethod::TAKE_SAME_HISTORY:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(*pairIndexes, nextSingDevMod,
                          resolveModule(window[i].module));
      if (localPair >= pairIndexes->size())
        continue;

      // Check if can be considered as the next CoincidenceEvent
      if (hist == window[i].hist) {
        result.iCoincidence = i;
        result.iPair = localPair;
        break;
      }
    }
    break;
  case CoincidenceMethod::TAKE_SAME_HISTORY_511:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(*pairIndexes, nextSingDevMod,
                          resolveModule(window[i].module));
      if (localPair >= pairIndexes->size())
        continue;

      // Check if can be considered as the next CoincidenceEvent
      double locald511 = std::fabs(window[i].e - 511.0);
      if (locald511 < d511 && hist == window[i].hist) {
        d511 = locald511;
        result.iCoincidence = i;
        result.iPair = localPair;
      }
    }
    break;
  case CoincidenceMethod::TAKE_CLOSEST:
  default:
    for (size_t i = 1; i < windowSize; ++i) {
      // Get pair
      localPair = getPair(*pairIndexes, nextSingDevMod,
                          resolveModule(window[i].module));
      if (localPair >= pairIndexes->size())
        continue;

      result.iCoincidence = i;
      result.iPair = localPair;
      break;
    }
    break;
  }

  return result;
}
