
#include "config/config.hh"
#include "common/constants.hh"
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <cstring>
#include <stdexcept>

SimConfig parseConfig(const std::string &filename) {

  std::ifstream infoS(filename);
  if (!infoS) {
    throw std::runtime_error("Unable to open information file '" + filename + "'");
  }

  SimConfig cfg;
  cfg.maxCounts = 10000000000000000ULL;
  cfg.header.version[0] = 7;
  cfg.header.version[1] = 1;
  cfg.header.calibrationID = 1;
  cfg.header.gatePeriod = 0;
  cfg.header.DOILayer = 1;
  cfg.header.method = 6;
  cfg.header.StudyID = 1;
  strcpy(cfg.header.identifier, "Simulation");
  cfg.header.startTime = 0.0;
  cfg.header.measurementTime = 0.0;
  cfg.header.weight = 1.0;
  cfg.header.maxTemp = 26.0;
  cfg.header.percentLoss = 0.0;

  std::string line;
  unsigned lineNum = 0;
  while (std::getline(infoS, line)) {
    lineNum++;
    size_t commentPos = line.find('#');
    if (commentPos != std::string::npos) {
      line = line.substr(0, commentPos);
    }
    
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    if (line.empty()) continue;

    size_t eqPos = line.find('=');
    if (eqPos == std::string::npos) {
      throw std::runtime_error("Invalid configuration line " + std::to_string(lineNum) + " (missing '='): " + line);
    }

    std::string key = line.substr(0, eqPos);
    std::string value = line.substr(eqPos + 1);
    
    key.erase(0, key.find_first_not_of(" \t"));
    key.erase(key.find_last_not_of(" \t") + 1);
    value.erase(0, value.find_first_not_of(" \t"));
    value.erase(value.find_last_not_of(" \t") + 1);

    if (key == "Counts") cfg.nCounts = std::stoull(value);
    else if (key == "DetectorSizeX") cfg.header.detectorSizeX = std::stod(value);
    else if (key == "DetectorSizeY") cfg.header.detectorSizeY = std::stod(value);
    else if (key == "DistanceBetweenDetectors") cfg.header.detectorDistance = std::stod(value);
    else if (key == "ModulesPerRing") cfg.header.moduleNumber = std::stoul(value);
    else if (key == "NumberOfRings") cfg.header.ringNumber = std::stoul(value);
    else if (key == "DistanceBetweenRings") cfg.header.ringDistance = std::stod(value);
    else if (key == "IsotopeName") strncpy(cfg.header.isotope, value.c_str(), sizeof(cfg.header.isotope) - 1);
    else if (key == "IsotopeHalfLife") cfg.header.isotopeHalfLife = std::stod(value);
    else if (key == "EnergyResolution") cfg.eRes = std::stod(value) / 100.0;
    else if (key == "PositionResolution") cfg.pRes = std::stod(value);
    else if (key == "TimeResolution") cfg.tRes = std::stod(value);
    else if (key == "CoincidenceWindow") cfg.coincidenceWindow = std::stod(value);
    else if (key == "EnergyWindow") cfg.energyWindow = std::stod(value) / 100.0;
    else if (key == "ModuleReconstructionBinsX") cfg.nBinsX = std::stoul(value);
    else if (key == "ModuleReconstructionBinsY") cfg.nBinsY = std::stoul(value);
    else if (key == "DetectorPixelsX") cfg.nBinsNormX = std::stoul(value);
    else if (key == "DetectorPixelsY") cfg.nBinsNormY = std::stoul(value);
    else if (key == "CrystalDepth") cfg.detectorDepth = std::stod(value);
    else if (key == "AcquisitionTime") cfg.header.acqTime = std::stod(value);
    else if (key == "OutputFormat") {
      int auxInt = std::stoi(value);
      if (auxInt == static_cast<int>(OutputFormat::LM_ONLY_COINCIDENCES)) {
        printf("Using LM format with only coincidences\n");
        cfg.outputFormat = OutputFormat::LM_ONLY_COINCIDENCES;
      } else {
        printf("Using LM format\n");
        cfg.outputFormat = OutputFormat::LM;
      }
    }
    else if (key == "CoincidenceMethod") {
      int auxInt = std::stoi(value);
      if (auxInt == static_cast<int>(CoincidenceMethod::TAKE_WINNER_IF_ALL_ARE_GOODS)) cfg.coinMethod = CoincidenceMethod::TAKE_WINNER_IF_ALL_ARE_GOODS;
      else if (auxInt == static_cast<int>(CoincidenceMethod::TAKE_CLOSEST)) cfg.coinMethod = CoincidenceMethod::TAKE_CLOSEST;
      else if (auxInt == static_cast<int>(CoincidenceMethod::TAKE_SAME_HISTORY)) cfg.coinMethod = CoincidenceMethod::TAKE_SAME_HISTORY;
      else if (auxInt == static_cast<int>(CoincidenceMethod::TAKE_SAME_HISTORY_511)) cfg.coinMethod = CoincidenceMethod::TAKE_SAME_HISTORY_511;
      else cfg.coinMethod = CoincidenceMethod::TAKE_CLOSEST;
    }
    else if (key == "ProjectionMethod") {
      int auxInt = std::stoi(value);
      if (auxInt == static_cast<int>(ProjectionMethod::EXTEND_DETECTOR)) cfg.projectionMethod = ProjectionMethod::EXTEND_DETECTOR;
      else if (auxInt == static_cast<int>(ProjectionMethod::FIT_IN_DETECTOR)) cfg.projectionMethod = ProjectionMethod::FIT_IN_DETECTOR;
      else cfg.projectionMethod = ProjectionMethod::NONE;
    }
    else if (key == "PairListFilename") cfg.pairListFilename = value;
    else if (key == "UseLogicalDetectors") cfg.useLogicalDetectors = (std::stoi(value) == 1);
    else if (key == "LogicalRings") cfg.nLogicalRings = std::stoul(value);
    else if (key == "DiscardMultiples") cfg.discardMultiples = std::stoi(value);
    else if (key == "GenerateHistogram") {
      int auxInt = std::stoi(value);
      if (auxInt == static_cast<int>(GenerateHistogram::LOR_INDEX)) cfg.generateHistogram = GenerateHistogram::LOR_INDEX;
      else cfg.generateHistogram = GenerateHistogram::NONE;
    }
    else if (key == "SaveWeight") cfg.saveWeight = (std::stoi(value) == 1);
    else if (key == "SaveMetadata") cfg.saveMetadata = (std::stoi(value) == 1);
    else {
      printf("Warning: Unknown configuration key '%s'\n", key.c_str());
    }
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
  
  cfg.header.rawCounts = static_cast<double>(cfg.nCounts);
  cfg.header.activity = activity * constants::BQ_TO_MICROCURIE;

  cfg.emin = 511.0e3 * (1.0 - cfg.energyWindow);
  cfg.emax = 511.0e3 * (1.0 + cfg.energyWindow);
  cfg.eminkeV = cfg.emin / 1.0e3;

  return cfg;
}
