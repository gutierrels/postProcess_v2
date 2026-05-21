
#ifndef POSTPROCESS_OUTPUT_LM_WRITER_HH
#define POSTPROCESS_OUTPUT_LM_WRITER_HH

#include "common/types.hh"
#include "config/config.hh"
#include <memory>
#include <vector>

class LMWriter {
public:
  LMWriter(const LMHeader &header, bool onlyCoincidences);
  ~LMWriter();

  void writeCoincidence(const CoincidenceEvent &c);

  void finalize();

  static std::unique_ptr<LMWriter> create(const SimConfig &cfg);

private:
  FILE *fLM;
  bool onlyCoincidences;
};



#endif
