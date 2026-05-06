
#ifndef POSTPROCESS_COINCIDENCE_ENGINE_HH
#define POSTPROCESS_COINCIDENCE_ENGINE_HH

#include "common/types.hh"
#include <vector>

struct CoincidenceResult {
  unsigned iCoincidence; // 0 = no coincidence
  unsigned iPair;
};

class CoincidenceEngine {
public:
  CoincidenceEngine(unsigned method, const std::vector<detPair> &pairs,
                    unsigned modPerRing, bool logicalMode);

  CoincidenceResult findCoincidence(const std::vector<single> &window,
                                    size_t windowSize) const;

private:
  unsigned coinMethod;
  const std::vector<detPair> &pairIndexes;
  unsigned modPerRing;
  bool useLogicalDetectors;
};

#endif
