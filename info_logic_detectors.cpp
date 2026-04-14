#include "common.hh"
#include <cstddef>
#include <cstdio>

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("usage: %s infoFile \n", argv[0]);
    return 1;
  }

  printf("LM Header size  : %lu\n",
         static_cast<unsigned long>(sizeof(LMHeader)));
  printf("Coincidence size: %lu\n",
         static_cast<unsigned long>(sizeof(coincidence)));

  std::ifstream infoS(argv[1]);

  if (!infoS) {
    printf("Unable to open information file '%s'\n", argv[1]);
    return -2;
  }

  LMHeader header;

  double eRes;
  double pRes; // cm

  double coincidenceWindow = 5.0e9; // s
  double energyWindow = 15.0;       // %
  unsigned nBinsX;
  unsigned nBinsY;

  unsigned nBinsNormX;
  unsigned nBinsNormY;

  double detectorDepth; // cm

  unsigned long long nCounts;
  unsigned long long maxCounts = 10000000000000000;

  unsigned outputFormat;
  unsigned projectionMethod;
  unsigned coinMethod;
  unsigned normMethod;
  unsigned attenuationMethod = ATTENUATION_METHOD::CYLINDER;

  std::string normFilename;
  std::string pairListFilename;

  double maxTimestamp = 0.0;

  // Read number of counts
  infoS >> nCounts;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read detector size in X axis (cm)
  infoS >> header.detectorSizeX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read detector size in Y axis (cm)
  infoS >> header.detectorSizeY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read distance between detectors (cm)
  infoS >> header.detectorDistance;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read number of modules per ring
  infoS >> header.moduleNumber;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read number of rings
  infoS >> header.ringNumber;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Read distance between rings (cm)
  infoS >> header.ringDistance;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Isotope name
  infoS >> header.isotope;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Isotope half life (s)
  infoS >> header.isotopeHalfLife;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Energy resolution (%)
  infoS >> eRes;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  eRes /= 100.0;

  // Position resolution (cm)
  infoS >> pRes;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Coincidence window (s)
  infoS >> coincidenceWindow;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Energy window (%)
  infoS >> energyWindow;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  energyWindow /= 100.0;

  // Module reconstruction bins in X axis
  infoS >> nBinsX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Module reconstruction bins in Y axis
  infoS >> nBinsY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Detector pixels in X axis
  infoS >> nBinsNormX;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Detector pixels in Y axis
  infoS >> nBinsNormY;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Crystal depth in cm
  infoS >> detectorDepth;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Acquisition time (s)
  infoS >> header.acqTime;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Save output format
  int auxInt;
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES:
    printf("Using Bruker LM format with only coincidences\n");
    outputFormat = OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES;
    break;
  case OUTPUT_FORMAT::CASTOR:
    printf("Using CASTOR LM format\n");
    outputFormat = OUTPUT_FORMAT::CASTOR;
    break;
  case OUTPUT_FORMAT::BRUKER_LM:
  default:
    printf("Using Bruker LM format\n");
    outputFormat = OUTPUT_FORMAT::BRUKER_LM;
    break;
  }

  // Coincidences method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS:
    coinMethod = COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS;
    break;
  case COINCIDENCE_METHOD::TAKE_CLOSEST:
    coinMethod = COINCIDENCE_METHOD::TAKE_CLOSEST;
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY:
    coinMethod = COINCIDENCE_METHOD::TAKE_SAME_HISTORY;
    break;
  case COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511:
    coinMethod = COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511;
    break;
  default:
    coinMethod = COINCIDENCE_METHOD::TAKE_CLOSEST;
  }

  // Projection method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case PROJECTION_METHOD::EXTEND_DETECTOR:
    projectionMethod = PROJECTION_METHOD::EXTEND_DETECTOR;
    break;
  case PROJECTION_METHOD::FIT_IN_DETECTOR:
    projectionMethod = PROJECTION_METHOD::FIT_IN_DETECTOR;
    break;
  default:
    projectionMethod = PROJECTION_METHOD::NONE;
  }

  // normalization method
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case NORMALIZATION_METHOD::CREATE:
    normMethod = NORMALIZATION_METHOD::CREATE;
    break;
  case NORMALIZATION_METHOD::APPLY:
    normMethod = NORMALIZATION_METHOD::APPLY;
    break;
  default:
    normMethod = NORMALIZATION_METHOD::NONE;
  }

  // normalization file
  infoS >> normFilename;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Attenuation correction
  infoS >> auxInt;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  switch (auxInt) {
  case ATTENUATION_METHOD::NONE:
    attenuationMethod = ATTENUATION_METHOD::NONE;
    break;
  case ATTENUATION_METHOD::CYLINDER:
    attenuationMethod = ATTENUATION_METHOD::CYLINDER;
    break;
  default:
    attenuationMethod = ATTENUATION_METHOD::NONE;
  }

  // Detector list filename

  infoS >> pairListFilename;
  infoS.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  if (!infoS) {
    printf("Corrupted information file '%s'\n", argv[1]);
    return -2;
  }
  infoS.close();

  double lambda = log(2) / header.isotopeHalfLife;
  double n0 = static_cast<double>(nCounts) / exp(-lambda * header.acqTime);
  double activity = lambda * n0;

  // Other LM Header parameters
  strcpy(header.identifier, "Simulation");
  header.rawCounts = static_cast<double>(nCounts);
  header.activity = activity * 2.7027027027027E-8; // In micro Curie
  header.startTime = 0.0;
  header.measurementTime = 0.0;
  header.weight = 1.0;
  header.maxTemp = 26.0;
  header.percentLoss = 0.0; //??
  // header.reserved = ;
  header.version[0] = 7;
  header.version[1] = 1;
  header.calibrationID = 1;
  header.gatePeriod = 0;
  header.DOILayer = 1;
  header.method = 6;
  header.StudyID = 1;

  // Physical detector
  const double physicalDetectorSizeX = header.detectorSizeX; // cm
  const double physicalDetectorSizeY = header.detectorSizeY; // cm

  const double physicalDBinSizeX =
      physicalDetectorSizeX / static_cast<double>(nBinsX); // cm

  const double physicalDBinSizeY =
      physicalDetectorSizeY / static_cast<double>(nBinsY); // cm

  const double physicalDBinSizeNormX =
      physicalDetectorSizeX / static_cast<double>(nBinsNormX);

  const double physicalDBinSizeNormY =
      physicalDetectorSizeY / static_cast<double>(nBinsNormY);

  const double physicalRingRadius = header.detectorDistance / 2.0;
  const double physicalRingDistance = header.ringDistance;

  const unsigned nPhysicalRings = header.ringNumber;
  const unsigned modPerRing = header.moduleNumber;
  const unsigned nPhysicalTotalMods = nPhysicalRings * modPerRing;

  printf("Physical detector size X: %f cm\n", physicalDetectorSizeX);
  printf("Physical detector size Y: %f cm\n", physicalDetectorSizeY);
  printf("Physical detector distance: %f cm\n", header.detectorDistance);
  printf("Physical detector distance between rings: %f cm\n",
         physicalRingDistance);
  printf("Physical detector rings: %u\n", nPhysicalRings);
  printf("Physical detector modules per ring: %u\n", modPerRing);
  printf("Physical detector total modules: %u\n", nPhysicalTotalMods);
  printf("Physical detector pixel size: %f x %f mm\n",
         physicalDBinSizeNormX * 10.0, physicalDBinSizeNormY * 10.0);
  printf("Physical detector bin size: %f x %f cm\n", physicalDBinSizeX,
         physicalDBinSizeY);
  printf("-------------------------------------\n");

  const double scannerLength = physicalDetectorSizeY * nPhysicalRings +
                               physicalRingDistance * nPhysicalRings;

  printf("Scanner length: %f cm (Contando medio gap arriba y medio abajo)\n",
         scannerLength);
  printf("-------------------------------------\n");

  // Logical detector
  const double logicalDetectorSizeY =
      scannerLength / (nPhysicalRings * 2);                 // cm
  const double logicalDetectorSizeX = logicalDetectorSizeY; // cm

  const double logicalDBinSizeX =
      logicalDetectorSizeX / static_cast<double>(nBinsX); // cm

  const double logicalDBinSizeY =
      logicalDetectorSizeY / static_cast<double>(nBinsY); // cm

  const double logicalDBinSizeNormX =
      logicalDetectorSizeX / static_cast<double>(nBinsNormX);

  const double logicalDBinSizeNormY =
      logicalDetectorSizeY / static_cast<double>(nBinsNormY);

  const double logicalRingRadius = physicalRingRadius; // cm
  const double logicalRingDistance = 0.0;

  const unsigned nLogicalRings = nPhysicalRings * 2;
  const unsigned nLogicalTotalMods =
      nLogicalRings * modPerRing; // modPerRing doesn't change

  printf("Logical detector size X: %f cm\n", logicalDetectorSizeX);
  printf("Logical detector size Y: %f cm\n", logicalDetectorSizeY);
  printf("Logical detector distance: %f cm\n", header.detectorDistance);
  printf("Logical detector distance between rings: %f cm\n",
         logicalRingDistance);
  printf("Logical detector rings: %u\n", nLogicalRings);
  printf("Logical detector modules per ring: %u\n", modPerRing);
  printf("Logical detector total modules: %u\n", nLogicalTotalMods);
  printf("Logical detector pixel size: %f x %f mm\n",
         logicalDBinSizeNormX * 10.0, logicalDBinSizeNormY * 10.0);
  printf("Logical detector bin size: %f x %f cm\n", logicalDBinSizeX,
         logicalDBinSizeY);

  printf("-------------------------------------\n");

  // Convert distances to mm
  header.detectorSizeX = logicalDetectorSizeX * 10.0;
  header.detectorSizeY = logicalDetectorSizeY * 10.0;
  header.detectorDistance *= 10.0;
  header.ringDistance = logicalRingDistance;
  header.ringNumber = nLogicalRings;

  header.print();

  return 0;
}
