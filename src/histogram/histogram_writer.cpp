
#include "histogram/histograms.hh"
#include <cstdio>
#include <filesystem>

void HistogramSet::writeAll(const std::string &outputDir, double maxTimestamp,
                            unsigned long nCoincidences) const {

  std::filesystem::create_directories(outputDir);


    std::string path = outputDir + "/histP.dat";
    FILE *fout = fopen(path.c_str(), "w");
    if (fout != nullptr) {
      for (unsigned v : histP)
        fprintf(fout, "%u\n", v);
      fclose(fout);
    }

    // Write Randoms per pair histogram
    path = outputDir + "/histR.dat";
    fout = fopen(path.c_str(), "w");
    if (fout != nullptr) {
      for (unsigned v : histR)
        fprintf(fout, "%u\n", v);
      fclose(fout);
    }

    for (size_t iPair = 0; iPair < nPairs; ++iPair) {

      // X

      // First pair module
      std::string histFilename =
          outputDir + "/histPx_" + std::to_string(iPair) + "_1.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < cfg.nBinsX; ++i)
          fprintf(fout, "%u\n", histPx[2 * iPair * cfg.nBinsX + i]);
        fclose(fout);
      }

      // Second pair module
      histFilename = outputDir + "/histPx_" + std::to_string(iPair) + "_2.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < cfg.nBinsX; ++i)
          fprintf(fout, "%u\n",
                  histPx[2 * iPair * cfg.nBinsX + cfg.nBinsX + i]);
        fclose(fout);
      }

      // Y

      // First pair module
      histFilename = outputDir + "/histPy_" + std::to_string(iPair) + "_1.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < cfg.nBinsY; ++i)
          fprintf(fout, "%u\n", histPy[2 * iPair * cfg.nBinsY + i]);
        fclose(fout);
      }

      // Second pair module
      histFilename = outputDir + "/histPy_" + std::to_string(iPair) + "_2.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < cfg.nBinsY; ++i)
          fprintf(fout, "%u\n",
                  histPy[2 * iPair * cfg.nBinsY + cfg.nBinsY + i]);
        fclose(fout);
      }
    }

    // Write LOR histogram
    if (cfg.generateHistogram == GENERATE_HISTOGRAM::LOR_INDEX) {
      fout = fopen("LOR.hist", "wb");
      if (fout != nullptr) {
        fwrite(histLOR.data(), sizeof(float), histLOR.size(), fout);
        fclose(fout);
      }
    }

    // Normalize to the highest central LOR value
    float histNormFactor = 1.0;
    std::vector<float> histNormCopy = histNorm;
    for (size_t i = 0; i < nPairs; ++i) {
      const size_t init = i * histNormPairSize;
      if (histNormFactor < histNormCopy[init + histNormCentralLOR])
        histNormFactor = histNormCopy[init + histNormCentralLOR];
    }
    for (float &e : histNormCopy) {
      e /= histNormFactor;
      if (e > 1000.0)
        e = 0.0;
    }
    printf("Normalization histogram factor: %E\n", histNormFactor);

    fout = fopen("measure.dat", "wb");
    if (fout != nullptr) {
      fwrite(histNormCopy.data(), sizeof(float), histNormCopy.size(), fout);
      fclose(fout);
    }


  writeModuleHistograms(outputDir);
}

void HistogramSet::writeModuleHistograms(const std::string &outputDir) const {
  // Write Module histograms
  std::string path = outputDir + "/histM.dat";
  FILE *fout = fopen(path.c_str(), "w");
  if (fout != nullptr) {
    for (unsigned v : histM)
      fprintf(fout, "%u\n", v);
    fclose(fout);
  }

  for (size_t imod = 0; imod < geo.totalMods; ++imod) {
    // X
    std::string histFilename =
        outputDir + "/histMx_" + std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < cfg.nBinsX; ++i)
        fprintf(fout, "%u\n", histMx[imod * cfg.nBinsX + i]);
      fclose(fout);
    }

    // Y
    histFilename = outputDir + "/histMy_" + std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < cfg.nBinsY; ++i)
        fprintf(fout, "%u\n", histMy[imod * cfg.nBinsY + i]);
      fclose(fout);
    }

    // E
    histFilename = outputDir + "/histMe_" + std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsE; ++i)
        fprintf(fout, "%E %u\n",
                geo.eminkeV + dekeV * (static_cast<double>(i) + 0.5),
                histMe[imod * nBinsE + i]);
      fclose(fout);
    }

    // E per z layer
    histFilename = outputDir + "/histMez_" + std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsE; ++i) {
        fprintf(fout, "%E ",
                geo.eminkeV + dekeV * (static_cast<double>(i) + 0.5));
        for (size_t j = 0; j < nBinsZ; ++j) {
          fprintf(fout, "%u ",
                  histMez[imod * nBinsE * nBinsZ + i * nBinsZ + j]);
        }
        fprintf(fout, "\n");
      }
      fclose(fout);
    }
  }
}
