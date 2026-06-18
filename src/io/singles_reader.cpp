
#include "io/singles_reader.hh"
#include "common/math_utils.hh"
#include "geometry/projection.hh"
#include <algorithm>

SinglesReader::SinglesReader(const std::string &prefix,
                             const DetectorGeometry &geo, const SimConfig &cfg)
    : geo(geo), cfg(cfg), done(false) {

  std::random_device rd{};
  gen = std::mt19937{rd()};
  normDist = std::normal_distribution<double>{0.0, 1.0};

  fmods.resize(geo.nPhysicalTotalMods);
  fileSizes.resize(geo.nPhysicalTotalMods);

  for (size_t i = 0; i < fmods.size(); ++i) {
    std::string dataFilename = prefix + std::to_string(i + 1) + ".dat";

    std::filesystem::path p = std::filesystem::current_path() / dataFilename;
    fileSizes[i] = std::filesystem::file_size(p);

    fmods[i] = fopen(dataFilename.c_str(), "rb");
    if (fmods[i] == nullptr) {
      throw std::runtime_error("Missing data file '" + dataFilename + "'");
    }
  }

  const double readingUpperLimit =
      (cfg.coinMethod == CoincidenceMethod::MULTIPLEXING)
          ? (cfg.promptGammaMax * constants::EV_PER_KEV)
          : cfg.emax;

  noInTimeSingles.resize(fmods.size());
  for (size_t i = 0; i < fmods.size(); ++i) {
    int err = noInTimeSingles[i].readPenRed(
        fmods[i], i, cfg.emin, readingUpperLimit, cfg.eRes, cfg.tRes,
        cfg.saveWeight, cfg.saveMetadata, normDist, gen);
    if (err != 0) {
      throw std::runtime_error("Error: Empty or corrupted module file (" +
                               std::to_string(i) + ")");
    }
  }

  if (cfg.useLogicalDetectors) {
    for (size_t i = 0; i < fmods.size(); i++) {
      toLogical(noInTimeSingles[i], geo.modPerRing, geo.logicalDetectorSizeY,
                geo.Rz, geo.physicalDz);
    }
  }
}

SinglesReader::~SinglesReader() {
  for (FILE *f : fmods) {
    if (f != nullptr) {
      fclose(f);
    }
  }
}

void SinglesReader::readNext(size_t fileIndex) {
  if (fmods[fileIndex] != nullptr) {
    const double readingUpperLimit =
        (cfg.coinMethod == CoincidenceMethod::MULTIPLEXING)
            ? (cfg.promptGammaMax * constants::EV_PER_KEV)
            : cfg.emax;
    int errRead = noInTimeSingles[fileIndex].readPenRed(
        fmods[fileIndex], fileIndex, cfg.emin, readingUpperLimit, cfg.eRes, cfg.tRes,
        cfg.saveWeight, cfg.saveMetadata, normDist, gen);

    if (errRead != 0) {
      if (errRead == 1) {
        // End of file
        fclose(fmods[fileIndex]);
        fmods[fileIndex] = nullptr;
      } else {
        throw std::runtime_error("Corrupted singles file for module " +
                                 std::to_string(fileIndex));
      }
    } else if (cfg.useLogicalDetectors) {
      toLogical(noInTimeSingles[fileIndex], geo.modPerRing,
                geo.logicalDetectorSizeY, geo.Rz, geo.physicalDz);
    }
  }
}

SingleEvent SinglesReader::seedFirst() {
  auto firstSingleIt =
      std::min_element(noInTimeSingles.begin(), noInTimeSingles.end());

  SingleEvent s = *firstSingleIt;

  unsigned fileIndex = cfg.useLogicalDetectors
                           ? static_cast<unsigned>(std::distance(
                                 noInTimeSingles.begin(), firstSingleIt))
                           : s.module;

  readNext(fileIndex);

  if (s.t >= constants::SENTINEL_TIME) {
    done = true;
  }
  return s;
}

void SinglesReader::fillWindow(std::vector<SingleEvent> &nextSingles,
                               double endTime, double timeMargin) {
  for (size_t i = 0; i < fmods.size(); ++i) {
    while (noInTimeSingles[i].true_t < endTime + timeMargin) {
      nextSingles.push_back(noInTimeSingles[i]);
      readNext(i);
    }
  }
}

bool SinglesReader::allDone() const { return done; }

std::string SinglesReader::progressString() const {
  std::string progressStr;
  for (size_t im = 0; im < fmods.size(); ++im) {
    char prefix[32];
    snprintf(prefix, sizeof(prefix), " Module %4zu", im + 1);
    progressStr += stringifyFileProgress(fmods[im], fileSizes[im], prefix);
  }
  return progressStr;
}
