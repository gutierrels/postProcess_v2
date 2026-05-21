
#include "histogram/histograms.hh"
#include <chrono>
#include <cstdio>
#include <filesystem>

void HistogramSet::writeAll(const std::string &outputDir,
                            const LMHeader &header) const {

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
        fprintf(fout, "%u\n", histPx[2 * iPair * cfg.nBinsX + cfg.nBinsX + i]);
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
        fprintf(fout, "%u\n", histPy[2 * iPair * cfg.nBinsY + cfg.nBinsY + i]);
      fclose(fout);
    }
  }

  // Write LOR histogram
  if (cfg.generateHistogram == GenerateHistogram::LOR_INDEX) {
    fout = fopen("LOR.hist", "wb");
    if (fout != nullptr) {
      // start timer
      auto start = std::chrono::high_resolution_clock::now();

      fwrite(&header, sizeof(LMHeader), 1, fout);
      fwrite(histLOR.data(), sizeof(float), histLOR.size(), fout);
      fclose(fout);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = end - start;
      printf("Time to write LOR histogram: %f s\n", elapsed.count());
    }
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
                cfg.eminkeV + dekeV * (static_cast<double>(i) + 0.5),
                histMe[imod * nBinsE + i]);
      fclose(fout);
    }

    // E per z layer
    histFilename = outputDir + "/histMez_" + std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsE; ++i) {
        fprintf(fout, "%E ",
                cfg.eminkeV + dekeV * (static_cast<double>(i) + 0.5));
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
