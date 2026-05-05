
#ifndef POSTPROCESS_OUTPUT_LM_WRITER_HH
#define POSTPROCESS_OUTPUT_LM_WRITER_HH

#include "common/types.hh"
#include "config/config.hh"
#include <memory>
#include <vector>

class LMWriter {
public:
  virtual ~LMWriter() = default;

  virtual void writeCoincidence(const coincidence &c) = 0;

  virtual void finalize() = 0;

  static std::unique_ptr<LMWriter> create(const SimConfig &cfg);
};

class BrukerLMWriter : public LMWriter {
public:
  BrukerLMWriter(const LMHeader &header, bool onlyCoincidences);
  ~BrukerLMWriter() override;

  void writeCoincidence(const coincidence &c) override;

  void finalize() override;

private:
  FILE *fLM;
  bool onlyCoincidences;
};



#endif
