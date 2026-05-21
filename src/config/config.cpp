
#include "config/config.hh"
#include "common/constants.hh"
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>

SimConfig parseConfig(const std::string &filename) {

  std::ifstream infoS(filename);
  if (!infoS) {
    throw std::runtime_error("Unable to open information file '" + filename + "'");
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
  case static_cast<int>(OutputFormat::BRUKER_LM_ONLY_COINCIDENCES):
    printf("Using Bruker LM format with only coincidences\n");
    cfg.outputFormat = OutputFormat::BRUKER_LM_ONLY_COINCIDENCES;
    break;

  case static_cast<int>(OutputFormat::BRUKER_LM):
  default:
    printf("Using Bruker LM format\n");
    cfg.outputFormat = OutputFormat::BRUKER_LM;
    break;
  }

  // Coincidences method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case static_cast<int>(CoincidenceMethod::TAKE_WINNER_IF_ALL_ARE_GOODS):
    cfg.coinMethod = CoincidenceMethod::TAKE_WINNER_IF_ALL_ARE_GOODS;
    break;
  case static_cast<int>(CoincidenceMethod::TAKE_CLOSEST):
    cfg.coinMethod = CoincidenceMethod::TAKE_CLOSEST;
    break;
  case static_cast<int>(CoincidenceMethod::TAKE_SAME_HISTORY):
    cfg.coinMethod = CoincidenceMethod::TAKE_SAME_HISTORY;
    break;
  case static_cast<int>(CoincidenceMethod::TAKE_SAME_HISTORY_511):
    cfg.coinMethod = CoincidenceMethod::TAKE_SAME_HISTORY_511;
    break;
  default:
    cfg.coinMethod = CoincidenceMethod::TAKE_CLOSEST;
  }

  // Projection method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case static_cast<int>(ProjectionMethod::EXTEND_DETECTOR):
    cfg.projectionMethod = ProjectionMethod::EXTEND_DETECTOR;
    break;
  case static_cast<int>(ProjectionMethod::FIT_IN_DETECTOR):
    cfg.projectionMethod = ProjectionMethod::FIT_IN_DETECTOR;
    break;
  default:
    cfg.projectionMethod = ProjectionMethod::NONE;
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
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case static_cast<int>(GenerateHistogram::LOR_INDEX):
    cfg.generateHistogram = GenerateHistogram::LOR_INDEX;
    break;
  default:
    cfg.generateHistogram = GenerateHistogram::NONE;
  }

  // PenRed Flags
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case 1:
    cfg.saveWeight = true;
    break;
  default:
    cfg.saveWeight = false;
  }

  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case 1:
    cfg.saveMetadata = true;
    break;
  default:
    cfg.saveMetadata = false;
  }

  if (!infoS) {
    throw std::runtime_error("Corrupted information file '" + filename + "'");
  }
  infoS.close();

  // Validate square detector pixels
  if (cfg.nBinsNormX == cfg.nBinsNormY) {
    cfg.N = cfg.nBinsNormX;
  } else {
    throw std::runtime_error("Not square detectors: '" + std::to_string(cfg.nBinsNormX) + "' and '" + std::to_string(cfg.nBinsNormY) + "' pixels");
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
  cfg.header.activity = activity * constants::BQ_TO_MICROCURIE;
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
