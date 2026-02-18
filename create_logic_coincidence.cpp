#include "common.hh"
#include <cstdio>

inline void toLogical(single &s, const unsigned modPerRing,
                      const double logicalDetectorSizeY) {}

int main(int argc, char **argv) {
  if (argc < 3) {
    printf("usage: %s infoFile input_prefix\n", argv[0]);
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
  header.detectorSizeX *= 10.0;
  header.detectorSizeY = logicalDetectorSizeY * 10.0;
  header.detectorDistance *= 10.0;
  header.ringDistance = logicalRingDistance;

  header.print();

  const double emin = 511.0e3 * (1.0 - energyWindow);
  const double emax = 511.0e3 * (1.0 + energyWindow);

  const double eminkeV = emin / 1.0e3;

  // Tallies
  unsigned long nCoincidences = 0;
  unsigned long nSingles = 0;

  // Create random gaussian distributions
  std::random_device rd{};
  std::mt19937 gen{rd()};

  std::normal_distribution normDist{0.0, 1.0};
  std::normal_distribution pBlur{0.0, pRes};

  // Create transformations for each module and ring
  std::vector<double> Dz(nLogicalRings);
  for (unsigned iRing = 0; iRing < nLogicalRings; ++iRing) {
    Dz[iRing] = iRing * (logicalRingDistance + logicalDetectorSizeY);
  }

  std::vector<std::array<double, 9>> Rz(modPerRing);
  constexpr double pi = 3.14159265359;
  const double angleStep = 2.0 * pi / static_cast<double>(modPerRing);

  for (unsigned im = 0; im < modPerRing; ++im) {

    // Calculate angle
    double angle = static_cast<double>(im) * angleStep;
    double cangle = cos(angle);
    double sangle = sin(angle);

    // Create rotation matrix
    const std::array<double, 9> localRz = {cangle, -sangle, 0.0, sangle, cangle,
                                           0.0,    0.0,     0.0, 1.0};

    Rz[im] = localRz;
  }

  // Read detector pair indexes
  std::vector<detPair> pairIndexes;
  FILE *fdetPair = fopen(pairListFilename.c_str(), "r");
  if (fdetPair == nullptr) {
    printf("Unable to open detector pair list file '%s'\n",
           pairListFilename.c_str());
    return -2;
  }

  char line[100];
  unsigned iline = 0;
  while (fgets(line, 100, fdetPair) != nullptr) {
    ++iline;
    unsigned a, b;
    if (sscanf(line, "%*u %u %u", &a, &b) != 2) {
      printf("Corrupted detector pair file '%s'. "
             "Check line %u: %s\n",
             pairListFilename.c_str(), iline, line);
      return -2;
    }
    pairIndexes.emplace_back(a, b);
  }
  fclose(fdetPair);
  printf("Loaded %lu detector pairs\n", pairIndexes.size());

  FILE *fLM;
  fLM = fopen("data.lm", "wb");
  if (fLM == nullptr) {
    printf("Unable to create LM file 'data.lm'\n");
    return -1;
  }

  // Write header
  fwrite(&header, sizeof(header), 1, fLM);

  // Open all module files
  std::vector<std::uintmax_t> fileSizes(nPhysicalTotalMods);
  std::vector<FILE *> fmods(nPhysicalTotalMods);
  for (size_t i = 0; i < fmods.size(); ++i) {

    std::string dataFilename(argv[2]);
    dataFilename += std::to_string(i + 1) + ".dat";

    // Get file size in bytes
    std::filesystem::path p = std::filesystem::current_path() / dataFilename;
    fileSizes[i] = std::filesystem::file_size(p);

    fmods[i] = fopen(dataFilename.c_str(), "rb");
    if (fmods[i] == nullptr) {
      printf("Missing data file '%s'\n", dataFilename.c_str());
      return -2;
    }
  }

  // Read first singles in data files
  std::vector<single> noInTimeSingles(fmods.size());
  for (size_t i = 0; i < fmods.size(); ++i) {
    int err = noInTimeSingles[i].readPenRed(fmods[i], i, emin, emax, eRes,
                                            normDist, gen);
    if (err != 0) {
      printf("Error: Empty or corrupted module file (%zu)\n", i);
      return -3;
    }
  }

  auto firstSingleIt =
      std::min_element(noInTimeSingles.begin(), noInTimeSingles.end());

  FILE *testFile = fopen("test_module.dat", "rb");

  single testSingle;
  testSingle.readPenRed(testFile, 0, emin, emax, eRes, normDist, gen);
  testSingle.stringify();
}