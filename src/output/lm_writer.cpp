
#include "output/lm_writer.hh"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

std::unique_ptr<LMWriter> LMWriter::create(const SimConfig &cfg) {
  bool onlyCoincidences =
      (cfg.outputFormat == OutputFormat::LM_ONLY_COINCIDENCES);
  return std::make_unique<LMWriter>(cfg.header, onlyCoincidences);
}

// ---------------------------------------------------------------------------
// LMWriter
// ---------------------------------------------------------------------------

LMWriter::LMWriter(const LMHeader &header, bool onlyCoincidences)
    : onlyCoincidences(onlyCoincidences) {
  fLM = fopen("data.lm", "wb");
  if (fLM == nullptr) {
    throw std::runtime_error("Unable to create LM file 'data.lm'");
  }
  // Write header
  fwrite(&header, sizeof(LMHeader), 1, fLM);
}

LMWriter::~LMWriter() {
  if (fLM != nullptr) {
    fclose(fLM);
    fLM = nullptr;
  }
}

void LMWriter::writeCoincidence(const CoincidenceEvent &c) {
  if (!onlyCoincidences || c.energy2 > 0.0) { // SingleEvent has energy2 = 0.0
    fwrite(&c, sizeof(CoincidenceEvent), 1, fLM);
  }
}

void LMWriter::finalize() {
  if (fLM != nullptr) {
    fclose(fLM);
    fLM = nullptr;
  }
}


