

#include "common.hh"

int main(int argc, char **argv) {

  if (argc < 3) {
    printf("usage: %s infoFile input_prefix\n", argv[0]);
    return 1;
  }

  printf("LM Header size  : %lu\n",
         static_cast<unsigned long>(sizeof(LMHeader)));
  printf("Coincidence size: %lu\n",
         static_cast<unsigned long>(sizeof(coincidence)));
  printf("Single size     : %lu\n", static_cast<unsigned long>(sizeof(single)));

  // ** Detector information ** //

  // Open and read information file
  std::ifstream infoS(argv[1]);
  if (!infoS) {
    printf("Unable to open information file '%s'\n", argv[1]);
    return -2;
  }

  // LM Header
  LMHeader header;

  double eRes;
  double pRes; // cm

  double coincidenceWindow; // s
  double energyWindow;      // %
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

  const double detectorSizeX = header.detectorSizeX;
  const double detectorSizeY = header.detectorSizeY;
  const double detectorSizeX05 = header.detectorSizeX / 2.0;
  const double detectorSizeY05 = header.detectorSizeY / 2.0;
  const double dBinSizeX = header.detectorSizeX / static_cast<double>(nBinsX);
  const double dBinSizeY = header.detectorSizeY / static_cast<double>(nBinsY);
  const double dBinSizeNormX =
      header.detectorSizeX / static_cast<double>(nBinsNormX);
  const double dBinSizeNormY =
      header.detectorSizeY / static_cast<double>(nBinsNormY);
  const double ringRad = header.detectorDistance / 2.0;
  const double ringDistance = header.ringDistance;

  // Convert distances to mm
  header.detectorSizeX *= 10.0;
  header.detectorSizeY *= 10.0;
  header.detectorDistance *= 10.0;
  header.ringDistance *= 10.0;

  const unsigned nRings = header.ringNumber;
  const unsigned modPerRing = header.moduleNumber;
  const unsigned totalMods = header.moduleNumber * header.ringNumber;

  //****************************//

  const double emin = 511.0e3 * (1.0 - energyWindow);
  const double emax = 511.0e3 * (1.0 + energyWindow);

  const double eminkeV = emin / 1.0e3;

  // Tallies
  unsigned long nCoincidences = 0;
  unsigned long nSingles = 0;
  unsigned long nRandoms = 0;
  unsigned long nTriples = 0;
  unsigned long nQuadruples = 0;
  unsigned long nQuintuplesOrMore = 0;

  // Create random gaussian distributions
  std::random_device rd{};
  std::mt19937 gen{rd()};

  std::normal_distribution normDist{0.0, 1.0};
  std::normal_distribution pBlur{0.0, pRes};

  // + Create transformations for each module and ring
  std::vector<double> Dz(nRings);
  for (unsigned iRing = 0; iRing < nRings; ++iRing) {
    Dz[iRing] = iRing * (ringDistance + detectorSizeY);
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

  header.print();

  // Read detector pair indexes
  std::vector<detPair> pairIndexes;
  FILE *fdetPair = fopen(pairListFilename.c_str(), "r");
  if (fdetPair == nullptr) {
    printf("Unable to open detector pair list file '%s'\n",
           pairListFilename.c_str());
    return -2;
  }

  printf("Reading detector pair list file '%s'\n", pairListFilename.c_str());

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

  // Create the LM file to be filled
  FILE *fLM;
  if (outputFormat != OUTPUT_FORMAT::CASTOR) {
    fLM = fopen("data.lm", "wb");
    if (fLM == nullptr) {
      printf("Unable to create LM file 'data.lm'\n");
      return -1;
    }
    // Write header
    fwrite(&header, sizeof(LMHeader), 1, fLM);
  }

  // Create Castor LM file to be filled
  FILE *fcLM;
  if (outputFormat == OUTPUT_FORMAT::CASTOR) {
    fcLM = fopen("cdata.lm", "wb");
    if (fcLM == nullptr) {
      printf("Unable to create LM file 'cdata.lm'\n");
      return -1;
    }
  }

  // Open all module files
  std::vector<std::uintmax_t> fileSizes(totalMods);
  std::vector<FILE *> fmods(totalMods);
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
      printf("Error: Empty or corrupted module file (%d)\n", i);
      return -3;
    }
  }

  // Get the first single
  auto firstSingleIt =
      std::min_element(noInTimeSingles.begin(), noInTimeSingles.end());

  // Save the first single in the processing vector
  std::vector<single> nextSingles;
  nextSingles.push_back(*firstSingleIt);

  // Read the next single to replace it
  unsigned firstSingleMod = firstSingleIt->module;
  int errRead = firstSingleIt->readPenRed(fmods[firstSingleMod], firstSingleMod,
                                          emin, emax, eRes, normDist, gen);
  if (errRead != 0) {
    if (errRead == 1) {
      // End of file
      fclose(fmods[firstSingleMod]);
      fmods[firstSingleMod] = nullptr;
    } else {
      // Corrupted file
      return -4;
    }
  }

  // Reserve histogram memory
  std::vector<unsigned> histP(pairIndexes.size(), 0);
  std::vector<unsigned> histR(pairIndexes.size(), 0);
  std::vector<unsigned> histPx(pairIndexes.size() * nBinsX * 2, 0);
  std::vector<unsigned> histPy(pairIndexes.size() * nBinsY * 2, 0);
  std::vector<unsigned> histM(totalMods, 0);
  std::vector<unsigned> histMx(totalMods * nBinsX, 0);
  std::vector<unsigned> histMy(totalMods * nBinsY, 0);
  const size_t nBinsE = 300;
  const double de = (emax - emin) / static_cast<double>(nBinsE);
  const double dekeV = de / 1.0e3;
  std::vector<unsigned> histMe(totalMods * nBinsE, 0);
  const size_t nBinsZ = 10;
  const double dz = detectorDepth / static_cast<double>(nBinsZ);
  std::vector<unsigned> histMez(nBinsZ * totalMods * nBinsE, 0);

  const size_t histNormPairSize =
      nBinsNormX * nBinsNormY * nBinsNormX * nBinsNormY;
  const size_t histNormCentralLOR =
      (nBinsNormX * nBinsNormY + 1) *
      (nBinsNormX * (nBinsNormY / 2) + nBinsNormX / 2);
  size_t histNormSize = histNormPairSize * pairIndexes.size();
  std::vector<float> histNorm(histNormSize, 0.0);

  std::map<std::pair<uint32_t, uint32_t>, double> cHistNorm;
  unsigned cHistNormSum = 0;
  double maxCountsLOR = 1.0;

  if (normMethod == NORMALIZATION_METHOD::APPLY) {
    // The normalization is applied, not created. Read factors for each LOR
    FILE *fnorm = fopen(normFilename.c_str(), "rb");
    if (fnorm == nullptr) {
      printf("Unable to read castor format normalization file '%s'\n",
             normFilename.c_str());
      return -3;
    }

    float factor;
    while (fread(&factor, sizeof(float), 1, fnorm) == 1) {
      uint32_t dets[2];
      if (fread(dets, sizeof(uint32_t), 2, fnorm) != 2) {
        printf("Corrupted castor format normalization file '%s'\n",
               normFilename.c_str());
        return -3;
      }
      std::pair<uint32_t, uint32_t> pair{dets[0], dets[1]};
      cHistNorm[pair] = factor;
    }
  }

  // Iterate over all singles
  size_t iterations = 0;
  std::string warnings;
  while (true) {

    // Get the initial window time
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
        }
      }
    }

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

    // Count multiples in the coincidence window
    if (nextSingles.size() == 3) {
      ++nTriples;
    } else if (nextSingles.size() == 4) {
      ++nQuadruples;
    } else if (nextSingles.size() >= 5) {
      ++nQuintuplesOrMore;
    }

    // Get first single history number
    unsigned long long hist = nextSingles[0].hist;

    // Find the closest single scored by a related detector
    double emaxSingle = 0.0;
    double d511 = 1.0e35;
    unsigned iCoincidence = 0;
    unsigned iPair = pairIndexes.size();
    unsigned nextSingDevMod = sim2devModule(modPerRing, nextSingles[0].module);
    unsigned localPair;
    switch (coinMethod) {
    case COINCIDENCE_METHOD::TAKE_WINNER_IF_ALL_ARE_GOODS:
      for (size_t i = 1; i < nextSingles.size(); ++i) {
        // Get pair
        localPair = getPair(pairIndexes, nextSingDevMod,
                            sim2devModule(modPerRing, nextSingles[i].module));
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
        // Get pair
        localPair = getPair(pairIndexes, nextSingDevMod,
                            sim2devModule(modPerRing, nextSingles[i].module));
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
        // Get pair
        localPair = getPair(pairIndexes, nextSingDevMod,
                            sim2devModule(modPerRing, nextSingles[i].module));
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
        // Get pair
        localPair = getPair(pairIndexes, nextSingDevMod,
                            sim2devModule(modPerRing, nextSingles[i].module));
        if (localPair >= pairIndexes.size())
          continue;

        // printf("%lu: %llu =? %llu\n ",
        //	 nextSingles.size(), hist, nextSingles[i].hist);
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
        p1 = toLocal(p1Aux, detectorSizeX, detectorSizeY, detectorDepth,
                     ringRad, Rz[nextSingles[0].module % modPerRing],
                     Dz[nextSingles[0].module / modPerRing], pBlur, gen);

        double p2Aux[3] = {nextSingles[iCoincidence].x,
                           nextSingles[iCoincidence].y,
                           nextSingles[iCoincidence].z};
        p2 = toLocal(p2Aux, detectorSizeX, detectorSizeY, detectorDepth,
                     ringRad, Rz[nextSingles[iCoincidence].module % modPerRing],
                     Dz[nextSingles[iCoincidence].module / modPerRing], pBlur,
                     gen);
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

        p1 =
            project(p1Blur, p2Blur, detectorSizeX, detectorSizeY, detectorDepth,
                    ringRad, Rz[nextSingles[0].module % modPerRing],
                    Dz[nextSingles[0].module / modPerRing],
                    projectionMethod == PROJECTION_METHOD::EXTEND_DETECTOR);

        p2 =
            project(p2Blur, p1Blur, detectorSizeX, detectorSizeY, detectorDepth,
                    ringRad, Rz[nextSingles[iCoincidence].module % modPerRing],
                    Dz[nextSingles[iCoincidence].module / modPerRing],
                    projectionMethod == PROJECTION_METHOD::EXTEND_DETECTOR);
      }

      // Check bounds
      if (p1[0] >= detectorSizeX || p1[0] <= 0.0 || p1[1] >= detectorSizeY ||
          p1[1] <= 0.0 || p1[2] >= detectorDepth || p1[2] < -1.0e-3 ||
          p2[0] >= detectorSizeX || p2[0] <= 0.0 || p2[1] >= detectorSizeY ||
          p2[1] <= 0.0 || p2[2] >= detectorDepth || p2[2] < -1.0e-3) {

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

        // Check if it is a random coincidence
        if (nextSingles[0].hist != nextSingles[iCoincidence].hist) {
          ++nRandoms;
          ++histR[iPair];
        }

        maxTimestamp = std::max(maxTimestamp, nextSingles[0].t);
        maxTimestamp = std::max(maxTimestamp, nextSingles[iCoincidence].t);

        if (outputFormat == OUTPUT_FORMAT::CASTOR) {
          // Castor format

          coincidenceCastor ccoin;
          std::pair<uint32_t, uint32_t> lor = castorLOR(
              nextSingles[0].module, nextSingles[iCoincidence].module, p1, p2,
              modPerRing, nBinsNormX, nBinsNormY, dBinSizeNormX, dBinSizeNormY);

          // Correct order (Consider symetry (d1,d2) = (d2,d1))
          if (pairIndexes[iPair].a !=
              sim2devModule(modPerRing, nextSingles[0].module)) {
            ccoin.det1 = lor.second;
            ccoin.det2 = lor.first;
          } else {
            ccoin.det1 = lor.first;
            ccoin.det2 = lor.second;
          }

          // Get timestamp
          ccoin.timestamp = static_cast<uint64_t>(nextSingles[0].t * 1.0e3);

          // Write timestamp
          fwrite(&(ccoin.timestamp), sizeof(uint32_t), 1, fcLM);

          if (attenuationMethod == ATTENUATION_METHOD::CYLINDER) {
            float attenuationFactor = attenuationFactorWithAir(
                lor0, lor1, 3.29841E-02, 2.96526E-02 * 1.20479000E-03,
                {0.0, 0.0, 7.5}, 3.5, 5.0);
            fwrite(&attenuationFactor, sizeof(float), 1, fcLM);
          }

          if (normMethod == NORMALIZATION_METHOD::CREATE) {
            auto it = cHistNorm.find(lor);
            if (it != cHistNorm.end()) {
              it->second += 1.0;
              if (maxCountsLOR < it->second)
                maxCountsLOR = it->second;
            } else {
              cHistNorm[lor] = 1.0;
            }
            ++cHistNormSum;
          } else if (normMethod == NORMALIZATION_METHOD::APPLY) {
            float factor = 0.0;
            auto it = cHistNorm.find(lor);
            if (it != cHistNorm.end()) {
              factor = it->second;
            }
            fwrite(&factor, sizeof(float), 1, fcLM);
          }

          fwrite(&(ccoin.det1), sizeof(uint32_t), 1, fcLM);
          fwrite(&(ccoin.det2), sizeof(uint32_t), 1, fcLM);
        } else {
          // Bruker based format

          // Check which pulse go first
          unsigned normBin1X, normBin1Y, normBin2X, normBin2Y;
          if (pairIndexes[iPair].a ==
              sim2devModule(modPerRing, nextSingles[0].module)) {
            c.energy1 = nextSingles[0].e;
            c.xPosition1 = static_cast<unsigned short>(p1[0] / dBinSizeX);
            c.yPosition1 = static_cast<unsigned short>(p1[1] / dBinSizeY);

            c.energy2 = nextSingles[iCoincidence].e;
            c.xPosition2 = static_cast<unsigned short>(p2[0] / dBinSizeX);
            c.yPosition2 = static_cast<unsigned short>(p2[1] / dBinSizeY);

            normBin1X = static_cast<unsigned>(p1[0] / dBinSizeNormX);
            normBin1Y = static_cast<unsigned>(p1[1] / dBinSizeNormY);

            normBin2X = static_cast<unsigned>(p2[0] / dBinSizeNormX);
            normBin2Y = static_cast<unsigned>(p2[1] / dBinSizeNormY);
          } else {
            c.energy1 = nextSingles[iCoincidence].e;
            c.xPosition1 = static_cast<unsigned short>(p2[0] / dBinSizeX);
            c.yPosition1 = static_cast<unsigned short>(p2[1] / dBinSizeY);

            c.energy2 = nextSingles[0].e;
            c.xPosition2 = static_cast<unsigned short>(p1[0] / dBinSizeX);
            c.yPosition2 = static_cast<unsigned short>(p1[1] / dBinSizeY);

            normBin1X = static_cast<unsigned>(p2[0] / dBinSizeNormX);
            normBin1Y = static_cast<unsigned>(p2[1] / dBinSizeNormY);

            normBin2X = static_cast<unsigned>(p1[0] / dBinSizeNormX);
            normBin2Y = static_cast<unsigned>(p1[1] / dBinSizeNormY);
          }

          if (normBin1X >= nBinsNormX) {
            printf("Unexpected bin index: %u - %u\n", nBinsNormX, normBin1X);
            fflush(stdout);
            return 1;
          }
          if (normBin2X >= nBinsNormX) {
            printf("Unexpected bin index: %u - %u\n", nBinsNormX, normBin2X);
            fflush(stdout);
            return 1;
          }

          if (normBin1Y >= nBinsNormY) {
            printf("Unexpected bin index: %u - %u\n", nBinsNormY, normBin1Y);
            fflush(stdout);
            return 1;
          }
          if (normBin2Y >= nBinsNormY) {
            printf("Unexpected bin index: %u - %u\n", nBinsNormY, normBin2Y);
            fflush(stdout);
            return 1;
          }

          normBin1X = nBinsNormX - 1 - normBin1X;
          normBin2X = nBinsNormX - 1 - normBin2X;

          c.pair = iPair;
          ++histP[iPair];
          ++histPx[2 * iPair * nBinsX + c.xPosition1];
          ++histPx[2 * iPair * nBinsX + nBinsX + c.xPosition2];
          ++histPy[2 * iPair * nBinsY + c.yPosition1];
          ++histPy[2 * iPair * nBinsY + nBinsY + c.yPosition2];

          const unsigned histPairIndex = iPair;
          const unsigned pairFirstIndex = histPairIndex * histNormPairSize;
          const unsigned binFirstIndex =
              (normBin1Y * nBinsNormX + normBin1X) * nBinsNormX * nBinsNormY;
          histNorm[pairFirstIndex + binFirstIndex + normBin2Y * nBinsNormX +
                   normBin2X] += 1.0;
        }

        ++histM[nextSingles[0].module];
        ++histM[nextSingles[iCoincidence].module];

        ++histMx[nextSingles[0].module * nBinsX + c.xPosition1];
        ++histMx[nextSingles[iCoincidence].module * nBinsX + c.xPosition2];

        ++histMy[nextSingles[0].module * nBinsY + c.yPosition1];
        ++histMy[nextSingles[iCoincidence].module * nBinsY + c.yPosition2];

        int ebin = (nextSingles[0].e - eminkeV) / dekeV;
        if (ebin >= 0 && ebin < nBinsE) {
          ++histMe[nextSingles[0].module * nBinsE + ebin];

          int zbin = static_cast<int>(p1[2] / dz);
          if (zbin >= 0 && zbin < nBinsZ) {
            ++histMez[nextSingles[0].module * nBinsE * nBinsZ + ebin * nBinsZ +
                      zbin];
          }
        }

        ebin = (nextSingles[iCoincidence].e - eminkeV) / dekeV;
        if (ebin >= 0 && ebin < nBinsE) {
          ++histMe[nextSingles[iCoincidence].module * nBinsE + ebin];

          int zbin = static_cast<int>(p2[2] / dz);
          if (zbin >= 0 && zbin < nBinsZ) {
            ++histMez[nextSingles[iCoincidence].module * nBinsE * nBinsZ +
                      ebin * nBinsZ + zbin];
          }
        }
      }

      // Erase used singles
      nextSingles.erase(nextSingles.begin() + iCoincidence);
      nextSingles.erase(nextSingles.begin());

    } else {
      // No coincidence detected for this single
      ++nSingles;

      // Get the point
      double p[3] = {nextSingles[0].x, nextSingles[0].y, nextSingles[0].z};

      // Apply the rotate to locate the point within the first ring detector
      matmul3D(Rz[nextSingles[0].module % modPerRing].data(), p);

      // Move the point to the first ring
      p[2] -= Dz[nextSingles[0].module / modPerRing];

      // Transform the point to local detector coordinates
      p[0] -= ringRad;
      p[1] += detectorSizeX05;
      p[2] = detectorSizeY - p[2];

      // Check bounds
      if (p[1] > detectorSizeX * 1.01 || p[1] < -1.0e-3 ||
          p[2] > detectorSizeY * 1.01 || p[2] < -1.0e-3 ||
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
      } else {

        c.energy1 = nextSingles[0].e;
        c.xPosition1 = static_cast<unsigned>(p[1] / dBinSizeX);
        c.yPosition1 = static_cast<unsigned>(p[2] / dBinSizeY);

        int ebin = (nextSingles[0].e - eminkeV) / dekeV;
        if (ebin >= 0 && ebin < nBinsE) {
          ++histMe[nextSingles[0].module * nBinsE + ebin];

          int zbin = static_cast<int>(p[0] / dz);
          if (zbin >= 0 && zbin < nBinsZ) {
            ++histMez[nextSingles[0].module * nBinsE * nBinsZ + ebin * nBinsZ +
                      zbin];
          }
        }

        c.energy2 = 0.0; // It is a single
        c.xPosition2 = 150;
        c.yPosition2 = 150;
        c.pair = nextSingles[0].module | 0x100;
      }
      // Erase used singles
      nextSingles.erase(nextSingles.begin());
    }

    // Save bruker format
    if (outputFormat != OUTPUT_FORMAT::CASTOR) {
      if (!outputFormat == OUTPUT_FORMAT::BRUKER_LM_ONLY_COINCIDENCES ||
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
      firstSingleMod = firstSingleIt->module;
      if (fmods[firstSingleMod] != nullptr) {
        errRead =
            firstSingleIt->readPenRed(fmods[firstSingleMod], firstSingleMod,
                                      emin, emax, eRes, normDist, gen);
        if (errRead != 0) {
          if (errRead == 1) {
            // End of file
            fclose(fmods[firstSingleMod]);
            fmods[firstSingleMod] = nullptr;
          } else {
            // Corrupted file
            return -4;
          }
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
        snprintf(prefix, 100, " Module %4d", im + 1);
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
    snprintf(prefix, 100, " Module %4d", im + 1);
    progressSstring += stringifyFileProgress(fmods[im], fileSizes[im], prefix);
  }
  std::cout << progressSstring << warnings << "\n";

  printf("Required iterations to process all files: %lu\n", iterations);
  printf("Number of coincidences: %lu\n", nCoincidences);
  printf("Number of singles     : %lu\n", nSingles);
  printf("Number of randoms     : %lu\n", nRandoms);
  printf("Number of triples     : %lu\n", nTriples);
  printf("Number of quadruples  : %lu\n", nQuadruples);
  printf("Number of 5+ multiples: %lu\n", nQuintuplesOrMore);

  if (outputFormat != OUTPUT_FORMAT::CASTOR) {
    fclose(fLM);

    // Write Pair histograms
    FILE *fout = fopen("modules/histP.dat", "w");
    if (fout != nullptr) {
      for (unsigned v : histP)
        fprintf(fout, "%u\n", v);
      fclose(fout);
    }

    // Write Randoms per pair histogram
    fout = fopen("modules/histR.dat", "w");
    if (fout != nullptr) {
      for (unsigned v : histR)
        fprintf(fout, "%u\n", v);
      fclose(fout);
    }

    for (size_t iPair = 0; iPair < pairIndexes.size(); ++iPair) {

      // X

      // First pair module
      std::string histFilename("modules/histPx_");
      histFilename += std::to_string(iPair) + "_1.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < nBinsX; ++i)
          fprintf(fout, "%u\n", histPx[2 * iPair * nBinsX + i]);
        fclose(fout);
      }

      // Second pair module
      histFilename.assign("modules/histPx_");
      histFilename += std::to_string(iPair) + "_2.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < nBinsX; ++i)
          fprintf(fout, "%u\n", histPx[2 * iPair * nBinsX + nBinsX + i]);
        fclose(fout);
      }

      // Y

      // First pair module
      histFilename.assign("modules/histPy_");
      histFilename += std::to_string(iPair) + "_1.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < nBinsY; ++i)
          fprintf(fout, "%u\n", histPy[2 * iPair * nBinsY + i]);
        fclose(fout);
      }

      // Second pair module
      histFilename.assign("modules/histPy_");
      histFilename += std::to_string(iPair) + "_2.dat";
      fout = fopen(histFilename.c_str(), "w");
      if (fout != nullptr) {
        for (size_t i = 0; i < nBinsY; ++i)
          fprintf(fout, "%u\n", histPy[2 * iPair * nBinsY + nBinsY + i]);
        fclose(fout);
      }
    }

    // Normalize to the highest central LOR value
    float histNormFactor = 1.0;
    for (size_t i = 0; i < pairIndexes.size(); ++i) {
      const size_t init = i * histNormPairSize;
      if (histNormFactor < histNorm[init + histNormCentralLOR])
        histNormFactor = histNorm[init + histNormCentralLOR];
    }
    for (float &e : histNorm) {
      e /= histNormFactor;
      if (e > 1000.0)
        e = 0.0;
    }
    printf("Normalization histogram factor: %E\n", histNormFactor);

    fout = fopen("measure.dat", "wb");
    fwrite(histNorm.data(), sizeof(float), histNorm.size(), fout);
    fclose(fout);
  } else {

    fclose(fcLM);

    if (normMethod == NORMALIZATION_METHOD::CREATE) {
      FILE *fout = fopen("cnorm.norm", "wb");

      const double normFact = 1.0 / static_cast<double>(cHistNormSum);
      // Get maximum LOR counts factor to re-normalize
      const double maxFact = maxCountsLOR * normFact;
      for (const auto &lor : cHistNorm) {
        const float factor = normFact * lor.second / maxFact;
        fwrite(&factor, sizeof(float), 1, fout);
        const uint32_t aux[2] = {lor.first.first, lor.first.second};
        fwrite(aux, sizeof(uint32_t), 2, fout);
      }
      fclose(fout);

      fout = fopen("cnorm.head", "w");

      fprintf(fout,
              "Scanner name                  : scanner\n"
              "Data filename                 : cnorm.norm\n"
              "Data mode                     : normalization\n"
              "Data type                     : PET\n"
              "Start time (s)                : 0\n"
              "Duration (s)                  : %f\n"
              "Number of events              : %lu\n"
              "Normalization correction flag : 1",
              maxTimestamp, static_cast<unsigned long>(cHistNorm.size()));

      fclose(fout);
    }

    FILE *fcH;
    fcH = fopen("cdata.cdh", "w");
    if (fcLM == nullptr) {
      printf("Unable to create LM header file 'cdata.cdh'\n");
    } else {
      fprintf(fcH,
              "Scanner name                  : scanner\n"
              "Data filename                 : cdata.lm\n"
              "Data mode                     : list-mode\n"
              "Data type                     : PET\n"
              "Start time (s)                : 0\n"
              "Duration (s)                  : %f\n"
              "Number of events              : %lu\n"
              "Attenuation correction flag   : %u\n"
              "Normalization correction flag : %u",
              maxTimestamp, static_cast<unsigned long>(nCoincidences),
              attenuationMethod == ATTENUATION_METHOD::CYLINDER ? 1 : 0,
              normMethod == NORMALIZATION_METHOD::APPLY ? 1 : 0);
      fclose(fcH);
    }
  }

  // Write Module histograms
  FILE *fout = fopen("modules/histM.dat", "w");
  if (fout != nullptr) {
    for (unsigned v : histM)
      fprintf(fout, "%u\n", v);
    fclose(fout);
  }

  for (size_t imod = 0; imod < totalMods; ++imod) {

    // X
    std::string histFilename("modules/histMx_");
    histFilename += std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsX; ++i)
        fprintf(fout, "%u\n", histMx[imod * nBinsX + i]);
      fclose(fout);
    }

    // Y
    histFilename.assign("modules/histMy_");
    histFilename += std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsY; ++i)
        fprintf(fout, "%u\n", histMy[imod * nBinsY + i]);
      fclose(fout);
    }

    // E
    histFilename.assign("modules/histMe_");
    histFilename += std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsE; ++i)
        fprintf(fout, "%E %u\n",
                eminkeV + dekeV * (static_cast<double>(i) + 0.5),
                histMe[imod * nBinsE + i]);
      fclose(fout);
    }

    // E per z layer
    histFilename.assign("modules/histMez_");
    histFilename += std::to_string(imod) + ".dat";
    fout = fopen(histFilename.c_str(), "w");
    if (fout != nullptr) {
      for (size_t i = 0; i < nBinsE; ++i) {
        fprintf(fout, "%E ", eminkeV + dekeV * (static_cast<double>(i) + 0.5));
        for (size_t j = 0; j < nBinsZ; ++j) {
          fprintf(fout, "%u ",
                  histMez[imod * nBinsE * nBinsZ + i * nBinsZ + j]);
        }
        fprintf(fout, "\n");
      }
      fclose(fout);
    }
  }

  return 0;
}
