
#ifndef POSTPROCESS_IO_SINGLES_READER_HH
#define POSTPROCESS_IO_SINGLES_READER_HH

#include "common/types.hh"
#include "config/config.hh"
#include "geometry/geometry.hh"
#include <filesystem>
#include <random>
#include <string>
#include <vector>

class SinglesReader {
public:
  SinglesReader(const std::string &prefix, const DetectorGeometry &geo,
                const SimConfig &cfg);
  ~SinglesReader();

  // Returns the first single to seed the process
  single seedFirst();

  // Fills nextSingles with all events in [t, t + windowEndTime]
  void fillWindow(std::vector<single> &nextSingles, double endTime,
                  double timeMargin);

  bool allDone() const;

  std::string progressString() const;

private:
  void readNext(size_t fileIndex);

  std::vector<FILE *> fmods;
  std::vector<std::uintmax_t> fileSizes;
  std::vector<single> noInTimeSingles;

  const DetectorGeometry &geo;
  const SimConfig &cfg;

  std::mt19937 gen;
  std::normal_distribution<double> normDist;

  bool done;
};

#endif
