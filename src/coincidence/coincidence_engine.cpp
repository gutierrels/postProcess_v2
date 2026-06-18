
#include "coincidence/coincidence_engine.hh"
#include "common/constants.hh"
#include "common/types.hh"
#include <cmath>

CoincidenceEngine::CoincidenceEngine(CoincidenceMethod method,
                                     const std::vector<DetectorPair> *pairs,
                                     unsigned modPerRing, bool logicalMode,
                                     double emin511KeV, double emax511KeV,
                                     double promptGammaMin,
                                     double promptGammaMax)
    : coinMethod(method), pairIndexes(pairs), modPerRing(modPerRing),
      useLogicalDetectors(logicalMode), emin511KeV(emin511KeV),
      emax511KeV(emax511KeV), promptGammaMin(promptGammaMin),
      promptGammaMax(promptGammaMax) {}

CoincidenceResult
CoincidenceEngine::findCoincidence(const std::vector<SingleEvent> &window,
                                   size_t windowSize) const {
  CoincidenceResult result = {0, static_cast<unsigned>(pairIndexes->size()),
                              false};

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

  case CoincidenceMethod::MULTIPLEXING:
    // window[0] is guaranteed to be a 511 keV photon (swap applied upstream),
    // but we double-check here to guard against any edge case.
    if (window[0].e >= emin511KeV && window[0].e <= emax511KeV) {
      for (size_t i = 1; i < windowSize; ++i) {
        // Only consider candidates that are also within the 511 keV window.
        // This prevents a prompt gamma from forming a false coincidence pair.
        if (window[i].e < emin511KeV || window[i].e > emax511KeV)
          continue;

        localPair = getPair(*pairIndexes, nextSingDevMod,
                            resolveModule(window[i].module));
        if (localPair >= pairIndexes->size())
          continue;

        // Found a valid 511-511 pair (closest in time)
        result.iCoincidence = i;
        result.iPair = localPair;

        // Scan remaining events for a prompt gamma
        for (size_t k = 0; k < windowSize; ++k) {
          if (k == 0 || k == i)
            continue;
          if (window[k].e >= promptGammaMin && window[k].e <= promptGammaMax) {
            result.isTriple = true;
            break;
          }
        }
        break;
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
