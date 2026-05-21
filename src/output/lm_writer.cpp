
#include "output/lm_writer.hh"
#include <cstdio>
#include <cstdlib>

std::unique_ptr<LMWriter> LMWriter::create(const SimConfig &cfg) {
  bool onlyCoincidences =
      (cfg.outputFormat == OutputFormat::BRUKER_LM_ONLY_COINCIDENCES);
  return std::make_unique<BrukerLMWriter>(cfg.header, onlyCoincidences);
}

// ---------------------------------------------------------------------------
// BrukerLMWriter
// ---------------------------------------------------------------------------

BrukerLMWriter::BrukerLMWriter(const LMHeader &header, bool onlyCoincidences)
    : onlyCoincidences(onlyCoincidences) {
  fLM = fopen("data.lm", "wb");
  if (fLM == nullptr) {
    printf("Unable to create LM file 'data.lm'\n");
    std::exit(-1);
  }
  // Write header
  fwrite(&header, sizeof(LMHeader), 1, fLM);
}

BrukerLMWriter::~BrukerLMWriter() {
  if (fLM != nullptr) {
    fclose(fLM);
    fLM = nullptr;
  }
}

void BrukerLMWriter::writeCoincidence(const CoincidenceEvent &c) {
  if (!onlyCoincidences || c.energy2 > 0.0) { // SingleEvent has energy2 = 0.0
    fwrite(&c, sizeof(CoincidenceEvent), 1, fLM);
  }
}



void BrukerLMWriter::finalize() {
  if (fLM != nullptr) {
    fclose(fLM);
    fLM = nullptr;
  }
}


