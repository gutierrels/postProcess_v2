
#include "config/config.hh"
#include <cstdio>
#include <fstream>
#include <limits>

SimConfig parseConfig(const std::string &filename) {

  std::ifstream infoS(filename);
  if (!infoS) {
    printf("Unable to open information file '%s'\n", filename.c_str());
    std::exit(-2);
  }

  SimConfig cfg;
  cfg.maxCounts = 10000000000000000ULL;

  // Read number of counts
  infoS >> cfg.nCounts;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read detector size in X axis (cm)
  infoS >> cfg.header.detectorSizeX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read detector size in Y axis (cm)
  infoS >> cfg.header.detectorSizeY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read distance between detectors (cm)
  infoS >> cfg.header.detectorDistance;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read number of modules per ring
  infoS >> cfg.header.moduleNumber;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read number of rings
  infoS >> cfg.header.ringNumber;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read distance between rings (cm)
  infoS >> cfg.header.ringDistance;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Isotope name
  infoS >> cfg.header.isotope;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Isotope half life (s)
  infoS >> cfg.header.isotopeHalfLife;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Energy resolution (%)
  infoS >> cfg.eRes;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  cfg.eRes /= 100.0;

  // Position resolution (cm)
  infoS >> cfg.pRes;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Time resolution (s)
  infoS >> cfg.tRes;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Coincidence window (s)
  infoS >> cfg.coincidenceWindow;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Energy window (%)
  infoS >> cfg.energyWindow;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  cfg.energyWindow /= 100.0;

  // Module reconstruction bins in X axis
  infoS >> cfg.nBinsX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Module reconstruction bins in Y axis
  infoS >> cfg.nBinsY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Detector pixels in X axis
  infoS >> cfg.nBinsNormX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Detector pixels in Y axis
  infoS >> cfg.nBinsNormY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Crystal depth in cm
  infoS >> cfg.detectorDepth;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Acquisition time (s)
  infoS >> cfg.header.acqTime;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Output format
  int auxInt;
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES:
    printf("Using Bruker LM format with only coincidences\n");
    cfg.outputFormat = OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES;
    break;

  case OUTPUT_FORMAT::BRUKER_LM:
  default:
    printf("Using Bruker LM format\n");
    cfg.outputFormat = OUTPUT_FORMAT::BRUKER_LM;
    break;
  }

  // Coincidences method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS:
    cfg.coinMethod = COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS;
    break;
  case COINCIDENCE_METHOD::TAKE_CLOSEST:
    cfg.coinMethod = COINCIDENCE_METHOD::TAKE_CLOSEST;
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY:
    cfg.coinMethod = COINCIDENCE_METHOD::TAKE_SAME_HISTORY;
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511:
    cfg.coinMethod = COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511;
    break;
  default:
    cfg.coinMethod = COINCIDENCE_METHOD::TAKE_CLOSEST;
  }

  // Projection method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case PROJECTION_METHOD::EXTEND_DETECTOR:
    cfg.projectionMethod = PROJECTION_METHOD::EXTEND_DETECTOR;
    break;
  case PROJECTION_METHOD::FIT_IN_DETECTOR:
    cfg.projectionMethod = PROJECTION_METHOD::FIT_IN_DETECTOR;
    break;
  default:
    cfg.projectionMethod = PROJECTION_METHOD::NONE;
  }


  // Detector list filename
  infoS >> cfg.pairListFilename;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Use logical detectors
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  cfg.useLogicalDetectors = (auxInt == 1);

  // Number of logical rings
  infoS >> cfg.nLogicalRings;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Discard multiples
  infoS >> cfg.discardMultiples;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Generate histogram
  infoS >> auxInt;
  switch (auxInt) {
  case GENERATE_HISTOGRAM::LOR_INDEX:
    cfg.generateHistogram = GENERATE_HISTOGRAM::LOR_INDEX;
    break;
  default:
    cfg.generateHistogram = GENERATE_HISTOGRAM::NONE;
  }

  if (!infoS) {
    printf("Corrupted information file '%s'\n", filename.c_str());
    std::exit(-2);
  }
  infoS.close();

  // Validate square detector pixels
  if (cfg.nBinsNormX == cfg.nBinsNormY) {
    cfg.N = cfg.nBinsNormX;
  } else {
    printf("Not square detectors: '%u' and '%u' pixels\n", cfg.nBinsNormX,
           cfg.nBinsNormY);
    std::exit(-2);
  }

  // Validate logical detector + CASTOR incompatibility
  if (cfg.useLogicalDetectors) {
    printf("Logical detectors enabled with %u logical rings\n",
           cfg.nLogicalRings);
  }

  // Populate remaining LMHeader fields
  double lambda = log(2) / cfg.header.isotopeHalfLife;
  double n0 =
      static_cast<double>(cfg.nCounts) / exp(-lambda * cfg.header.acqTime);
  double activity = lambda * n0;

  strcpy(cfg.header.identifier, "Simulation");
  cfg.header.rawCounts = static_cast<double>(cfg.nCounts);
  cfg.header.activity = activity * 2.7027027027027E-8; // In micro Curie
  cfg.header.startTime = 0.0;
  cfg.header.measurementTime = 0.0;
  cfg.header.weight = 1.0;
  cfg.header.maxTemp = 26.0;
  cfg.header.percentLoss = 0.0;
  cfg.header.version[0] = 7;
  cfg.header.version[1] = 1;
  cfg.header.calibrationID = 1;
  cfg.header.gatePeriod = 0;
  cfg.header.DOILayer = 1;
  cfg.header.method = 6;
  cfg.header.StudyID = 1;

  cfg.emin = 511.0e3 * (1.0 - cfg.energyWindow);
  cfg.emax = 511.0e3 * (1.0 + cfg.energyWindow);
  cfg.eminkeV = cfg.emin / 1.0e3;

  return cfg;
}
