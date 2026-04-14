#include "common.hh"
#include <cstddef>
#include <cstdio>

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
      scannerLength / (nPhysicalRings * 2);                 // entre 6     // cm
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
  std::vector<double> physicalDz(nPhysicalRings);
  for (unsigned iRing = 0; iRing < nPhysicalRings; ++iRing) {
    physicalDz[iRing] = iRing * (physicalRingDistance + physicalDetectorSizeY);
  }

  std::vector<double> logicalDz(nLogicalRings);
  for (unsigned iRing = 0; iRing < nLogicalRings; ++iRing) {
    logicalDz[iRing] = iRing * (logicalRingDistance + logicalDetectorSizeY);
  }

  // Rz is the same logical/physical

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

  for (size_t i = 0; i < fmods.size(); i++) {
    toLogical(noInTimeSingles[i], modPerRing, logicalDetectorSizeY, Rz,
              physicalDz);
  }

  auto firstSingleIt =
      std::min_element(noInTimeSingles.begin(), noInTimeSingles.end());

  // Save the first single in the processing vector
  std::vector<single> nextSingles;
  nextSingles.push_back(*firstSingleIt);

  // Read the next single to replace it
  unsigned firstSingleMod = firstSingleIt->module;

  unsigned fileIndex = std::distance(noInTimeSingles.begin(), firstSingleIt);

  int errRead = firstSingleIt->readPenRed(fmods[fileIndex], fileIndex, emin,
                                          emax, eRes, normDist, gen);
  if (errRead != 0) {
    if (errRead == 1) {
      // End of file
      fclose(fmods[fileIndex]);
      fmods[fileIndex] = nullptr;
    } else {
      // Corrupted file
      return -4;
    }
  } else {
    // Convert to logical detector
    toLogical(*firstSingleIt, modPerRing, logicalDetectorSizeY, Rz, physicalDz);
  }

  printf("-------------------------------------\n");

  // Reserve histogram memory

  std::vector<unsigned> histP(pairIndexes.size(), 0);
  std::vector<unsigned> histM(nLogicalTotalMods, 0);

  // Iterate over all singles
  size_t iterations = 0;
  std::string warnings;

  while (true) {

    double initTime = nextSingles[0].t;
    double endTime = initTime + coincidenceWindow;

    // Check if all singles have been processed
    if (initTime > 1.0e20)
      break; // End of execution

    // Read all singles within time window
    for (size_t i = 0; i < fmods.size(); ++i) {
      while (noInTimeSingles[i].t < endTime) {
        nextSingles.push_back(noInTimeSingles[i]);
        errRead = noInTimeSingles[i].readPenRed(fmods[i], i, emin, emax, eRes,
                                                normDist, gen);
        if (errRead != 0) {
          if (errRead == 1) {
            // End of file
            fclose(fmods[i]);
            fmods[i] = nullptr;
          } else {
            // Corrupted file
            return -4;
          }
        } else {
          // Convert to logical detector
          toLogical(noInTimeSingles[i], modPerRing, logicalDetectorSizeY, Rz,
                    physicalDz);
        }
      }
    }

    // Debug, delete later
    if (nextSingles.size() > 10000) {
      size_t nSing = nextSingles.size();
      printf("UPS start!\n"
             "%s\n%s\n%s\n%s\n%s\n"
             ".\n.\n.\n"
             "%s\n%s\n%s\n%s\n%s\n"
             "UPS end!\n",
             nextSingles[0].stringify().c_str(),
             nextSingles[1].stringify().c_str(),
             nextSingles[2].stringify().c_str(),
             nextSingles[3].stringify().c_str(),
             nextSingles[4].stringify().c_str(),
             nextSingles[nSing].stringify().c_str(),
             nextSingles[nSing - 1].stringify().c_str(),
             nextSingles[nSing - 2].stringify().c_str(),
             nextSingles[nSing - 3].stringify().c_str(),
             nextSingles[nSing - 4].stringify().c_str());
      return -1;
    }

    // Sort new singles
    std::sort(nextSingles.begin(), nextSingles.end());

    // Get first single history number
    unsigned long long hist = nextSingles[0].hist;

    // Find the closest single scored by a related logical detector
    double emaxSingle = 0.0;
    unsigned iCoincidence = 0;
    double d511 = 1.0e35;
    unsigned iPair = pairIndexes.size();
    // Use logical module index directly (already converted by toLogical)
    unsigned nextSingLogMod = nextSingles[0].module;
    unsigned localPair;

    switch (coinMethod) {

    case COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS:
      for (size_t i = 1; i < nextSingles.size(); ++i) {
        // Get pair using logical module indices directly
        localPair = getPair(pairIndexes, nextSingLogMod, nextSingles[i].module);
        if (localPair >= pairIndexes.size())
          continue;

        // Check if can be considered as the next coincidence
        if (emaxSingle < nextSingles[i].e) {
          emaxSingle = nextSingles[i].e;
          iCoincidence = i;
          iPair = localPair;
        }
      }
      break;

    case COINCIDENCE_METHOD::TAKE_SAME_HISTORY:
      for (size_t i = 1; i < nextSingles.size(); ++i) {
        // Get pair using logical module indices directly
        localPair = getPair(pairIndexes, nextSingLogMod, nextSingles[i].module);
        if (localPair >= pairIndexes.size())
          continue;

        // Check if can be considered as the next coincidence
        if (hist == nextSingles[i].hist) {
          iCoincidence = i;
          iPair = localPair;
          break;
        }
      }
      break;

    case COINCIDENCE_METHOD::TAKE_SAME_HISTORY_511:
      for (size_t i = 1; i < nextSingles.size(); ++i) {
        // Get pair using logical module indices directly
        localPair = getPair(pairIndexes, nextSingLogMod, nextSingles[i].module);
        if (localPair >= pairIndexes.size())
          continue;

        // Check if can be considered as the next coincidence
        double locald511 = std::fabs(nextSingles[i].e - 511.0);
        if (locald511 < d511 && hist == nextSingles[i].hist) {
          d511 = locald511;
          iCoincidence = i;
          iPair = localPair;
        }
      }
      break;

    case COINCIDENCE_METHOD::TAKE_CLOSEST:
    default:
      for (size_t i = 1; i < nextSingles.size(); ++i) {
        // Get pair using logical module indices directly
        localPair = getPair(pairIndexes, nextSingLogMod, nextSingles[i].module);
        if (localPair >= pairIndexes.size())
          continue;

        iCoincidence = i;
        iPair = localPair;
        break;
      }
      break;
    }

    // Create the entry
    coincidence c;
    c.time = nextSingles[0].t;
    c.amount = 1.0;
    c.gate_flag = 0;

    // Flag if a coincidence has been accepted
    bool acceptedCoincidence = false;

    // Check if some other single is in energy window to create a coincidence
    if (iCoincidence > 0) {

      // It is a possible coincidence

      // + Project both points to the detector surface + //

      // First point

      std::array<double, 3> p1;
      std::array<double, 3> p2;
      std::array<double, 3> lor0;
      std::array<double, 3> lor1;

      if (projectionMethod == PROJECTION_METHOD::NONE) {

        lor0 = std::array<double, 3>{nextSingles[0].x, nextSingles[0].y,
                                     nextSingles[0].z};
        lor1 = std::array<double, 3>{nextSingles[iCoincidence].x,
                                     nextSingles[iCoincidence].y,
                                     nextSingles[iCoincidence].z};

        double p1Aux[3] = {nextSingles[0].x, nextSingles[0].y,
                           nextSingles[0].z};
        p1 = toLocal(p1Aux, logicalDetectorSizeX, logicalDetectorSizeY,
                     detectorDepth, logicalRingRadius,
                     Rz[nextSingles[0].module % modPerRing],
                     logicalDz[nextSingles[0].module / modPerRing], pBlur, gen);

        double p2Aux[3] = {nextSingles[iCoincidence].x,
                           nextSingles[iCoincidence].y,
                           nextSingles[iCoincidence].z};
        p2 = toLocal(p2Aux, logicalDetectorSizeX, logicalDetectorSizeY,
                     detectorDepth, logicalRingRadius,
                     Rz[nextSingles[iCoincidence].module % modPerRing],
                     logicalDz[nextSingles[iCoincidence].module / modPerRing],
                     pBlur, gen);
      } else {

        // Apply bluring
        double p1Blur[3] = {
            nextSingles[0].x + pBlur(gen),
            nextSingles[0].y + pBlur(gen),
            nextSingles[0].z + pBlur(gen),
        };

        double p2Blur[3] = {
            nextSingles[iCoincidence].x + pBlur(gen),
            nextSingles[iCoincidence].y + pBlur(gen),
            nextSingles[iCoincidence].z + pBlur(gen),
        };

        lor0 = std::array<double, 3>{p1Blur[0], p1Blur[1], p1Blur[2]};
        lor1 = std::array<double, 3>{p2Blur[0], p2Blur[1], p2Blur[2]};

        p1 = project(p1Blur, p2Blur, logicalDetectorSizeX, logicalDetectorSizeY,
                     detectorDepth, logicalRingRadius,
                     Rz[nextSingles[0].module % modPerRing],
                     logicalDz[nextSingles[0].module / modPerRing],
                     projectionMethod == PROJECTION_METHOD::EXTEND_DETECTOR);

        p2 = project(p2Blur, p1Blur, logicalDetectorSizeX, logicalDetectorSizeY,
                     detectorDepth, logicalRingRadius,
                     Rz[nextSingles[iCoincidence].module % modPerRing],
                     logicalDz[nextSingles[iCoincidence].module / modPerRing],
                     projectionMethod == PROJECTION_METHOD::EXTEND_DETECTOR);
      }

      // Check bounds
      if (p1[0] >= logicalDetectorSizeX || p1[0] <= 0.0 ||
          p1[1] >= logicalDetectorSizeY || p1[1] <= 0.0 ||
          p1[2] >= detectorDepth || p1[2] < -1.0e-3 ||
          p2[0] >= logicalDetectorSizeX || p2[0] <= 0.0 ||
          p2[1] >= logicalDetectorSizeY || p2[1] <= 0.0 ||
          p2[2] >= detectorDepth || p2[2] < -1.0e-3) {

        if (projectionMethod != PROJECTION_METHOD::EXTEND_DETECTOR) {
          char auxBuffer[1000];
          snprintf(auxBuffer, 1000,
                   "Warning: Unconsistent coincidence data.\n"
                   " Original P1: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                   " Original P2: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                   "Resulting P1: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                   "Resulting P2: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                   "P1 Ring: %d, Module: %d (%d)\n"
                   "P2 Ring: %d, Module: %d (%d)\n",
                   nextSingles[0].e, nextSingles[0].x, nextSingles[0].y,
                   nextSingles[0].z, nextSingles[0].t,
                   nextSingles[iCoincidence].e, nextSingles[iCoincidence].x,
                   nextSingles[iCoincidence].y, nextSingles[iCoincidence].z,
                   nextSingles[iCoincidence].t, nextSingles[0].e, p1[0], p1[1],
                   p1[2], nextSingles[0].t, nextSingles[iCoincidence].e, p2[0],
                   p2[1], p2[2], nextSingles[iCoincidence].t,
                   nextSingles[0].module / modPerRing,
                   nextSingles[0].module % modPerRing, nextSingles[0].module,
                   nextSingles[iCoincidence].module / modPerRing,
                   nextSingles[iCoincidence].module % modPerRing,
                   nextSingles[iCoincidence].module);

          warnings += auxBuffer;
        }
      } else {

        // It is a confirmed coincidence
        acceptedCoincidence = true;
        ++nCoincidences;

        maxTimestamp = std::max(maxTimestamp, nextSingles[0].t);
        maxTimestamp = std::max(maxTimestamp, nextSingles[iCoincidence].t);

        c.pair = iPair;

        // Check which single corresponds to detector 'a' of the pair
        if (pairIndexes[iPair].a == nextSingles[0].module) {
          c.energy1 = nextSingles[0].e;
          c.xPosition1 = static_cast<unsigned short>(p1[0] / logicalDBinSizeX);
          c.yPosition1 = static_cast<unsigned short>(p1[1] / logicalDBinSizeY);

          c.energy2 = nextSingles[iCoincidence].e;
          c.xPosition2 = static_cast<unsigned short>(p2[0] / logicalDBinSizeX);
          c.yPosition2 = static_cast<unsigned short>(p2[1] / logicalDBinSizeY);
        } else {
          c.energy1 = nextSingles[iCoincidence].e;
          c.xPosition1 = static_cast<unsigned short>(p2[0] / logicalDBinSizeX);
          c.yPosition1 = static_cast<unsigned short>(p2[1] / logicalDBinSizeY);

          c.energy2 = nextSingles[0].e;
          c.xPosition2 = static_cast<unsigned short>(p1[0] / logicalDBinSizeX);
          c.yPosition2 = static_cast<unsigned short>(p1[1] / logicalDBinSizeY);
        }

        ++histP[iPair];

        ++histM[nextSingles[0].module];
        ++histM[nextSingles[iCoincidence].module];

        // Erase used singles
        nextSingles.erase(nextSingles.begin() + iCoincidence);
        nextSingles.erase(nextSingles.begin());
      }
    } else {

      // No coincidence detected for this single
      ++nSingles;

      // Get the point
      double p[3] = {nextSingles[0].x, nextSingles[0].y, nextSingles[0].z};

      // Apply the rotate to locate the point within the first ring detector
      matmul3D(Rz[nextSingles[0].module % modPerRing].data(), p);

      // Move the point to the first ring
      p[2] -= logicalDz[nextSingles[0].module / modPerRing];

      // Transform the point to local detector coordinates
      p[0] -= logicalRingRadius;
      p[1] += logicalDetectorSizeX / 2.0;
      p[2] = logicalDetectorSizeY - p[2];

      // Check bounds
      if (p[1] > logicalDetectorSizeX * 1.01 || p[1] < -1.0e-3 ||
          p[2] > logicalDetectorSizeY * 1.01 || p[2] < -1.0e-3 ||
          p[0] > detectorDepth * 1.01 || p[0] < -1.0e-3) {

        char auxBuffer[1000];
        snprintf(auxBuffer, 1000,
                 "Warning: Unconsistent single data.\n"
                 " Original P: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                 "Resulting P: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                 "Ring: %d, Module: %d (%d)\n",
                 nextSingles[0].e, nextSingles[0].x, nextSingles[0].y,
                 nextSingles[0].z, nextSingles[0].t, nextSingles[0].e, p[1],
                 p[2], p[0], nextSingles[0].t,
                 nextSingles[0].module / modPerRing,
                 nextSingles[0].module % modPerRing, nextSingles[0].module);
        warnings += auxBuffer;
      }

      // Erase used singles
      nextSingles.erase(nextSingles.begin());
    }

    if (outputFormat != OUTPUT_FORMAT::CASTOR) {
      if (outputFormat != OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES ||
          acceptedCoincidence)
        fwrite(&c, sizeof(coincidence), 1, fLM);
    }

    if (nSingles + nCoincidences >= maxCounts)
      break;

    // If the processing vector is empty, put in the next in time
    if (nextSingles.size() == 0) {

      // Get the first single
      firstSingleIt =
          std::min_element(noInTimeSingles.begin(), noInTimeSingles.end());

      // Save the first single in the processing vector
      nextSingles.push_back(*firstSingleIt);

      // Read the next single to replace it
      fileIndex = std::distance(noInTimeSingles.begin(), firstSingleIt);
      if (fmods[fileIndex] != nullptr) {
        errRead = firstSingleIt->readPenRed(fmods[fileIndex], fileIndex, emin,
                                            emax, eRes, normDist, gen);
        if (errRead != 0) {
          if (errRead == 1) {
            // End of file
            fclose(fmods[fileIndex]);
            fmods[fileIndex] = nullptr;
          } else {
            // Corrupted file
            return -4;
          }
        } else {
          // Convert to logical detector
          toLogical(*firstSingleIt, modPerRing, logicalDetectorSizeY, Rz,
                    physicalDz);
        }
      }
    }

    // Print progress
    if (iterations % 100000 == 0) {
      std::cout << "\033[2J";   // Clear window
      std::cout << "\033[1;1H"; // Go to first line
      std::string progressSstring;
      for (size_t im = 0; im < fmods.size(); ++im) {
        char prefix[20];
        snprintf(prefix, 20, " Module %4zu", im + 1);
        progressSstring +=
            stringifyFileProgress(fmods[im], fileSizes[im], prefix);
      }
      progressSstring += "Simultaneous processing singles: " +
                         std::to_string(nextSingles.size()) + "\n";

      std::cout << progressSstring << warnings << "\n";
    }
    ++iterations;
  }

  // Print final progress
  // std::cout << "\033[2J"; //Clear window
  // std::cout << "\033[1;1H"; //Go to first line
  std::string progressSstring;
  for (size_t im = 0; im < fmods.size(); ++im) {
    char prefix[20];
    snprintf(prefix, 20, " Module %4zu", im + 1);
    progressSstring += stringifyFileProgress(fmods[im], fileSizes[im], prefix);
  }
  std::cout << progressSstring << warnings << "\n";

  printf("Required iterations to process all files: %lu\n", iterations);
  printf("Number of coincidences: %lu\n", nCoincidences);
  printf("Number of singles     : %lu\n", nSingles);

  fclose(fLM);

  // Write Pair histograms
  FILE *fout = fopen("modules/histP.dat", "w");
  if (fout != nullptr) {
    for (unsigned v : histP)
      fprintf(fout, "%u\n", v);
    fclose(fout);
  }

  // Write Module histograms
  fout = fopen("modules/histM.dat", "w");
  if (fout != nullptr) {
    for (unsigned v : histM)
      fprintf(fout, "%u\n", v);
    fclose(fout);
  }

  header.print();

  return 0;
}