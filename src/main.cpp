
#include "coincidence/coincidence_engine.hh"
#include "common/constants.hh"
#include "common/math_utils.hh"
#include "common/types.hh"
#include "config/config.hh"
#include "geometry/geometry.hh"
#include "geometry/projection.hh"
#include "histogram/histograms.hh"
#include "io/pair_list.hh"
#include "io/singles_reader.hh"
#include "output/lm_writer.hh"


#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc < 3) {
    printf("usage: %s infoFile input_prefix\n", argv[0]);
    return 1;
  }

  printf("LM Header size  : %lu\n",
         static_cast<unsigned long>(sizeof(LMHeader)));
  printf("Coincidence size: %lu\n",
         static_cast<unsigned long>(sizeof(CoincidenceEvent)));
  printf("Single size     : %lu\n", static_cast<unsigned long>(sizeof(SingleEvent)));

  // Setup
  SimConfig cfg = parseConfig(argv[1]);
  cfg.header.print();

  DetectorGeometry geo = DetectorGeometry::build(cfg);
  std::vector<DetectorPair> pairs = readPairList(cfg.pairListFilename);
  SinglesReader reader(argv[2], geo, cfg);
  CoincidenceEngine engine(cfg.coinMethod, &pairs, geo.modPerRing,
                           cfg.useLogicalDetectors);
  HistogramSet hists(cfg, geo, pairs.size());
  std::unique_ptr<LMWriter> writer = LMWriter::create(cfg);

  // random bluring
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> pBlur{0.0, cfg.pRes / 2.355};

  // Tallies
  unsigned long nCoincidences = 0;
  unsigned long nSingles = 0;
  unsigned long nRandoms = 0;
  unsigned long nTriples = 0;
  unsigned long nQuadruples = 0;
  unsigned long nQuintuplesOrMore = 0;
  double maxTimestamp = 0.0;

  // Main loop
  std::vector<SingleEvent> window;
  window.push_back(reader.seedFirst());

  size_t iterations = 0;
  std::string warnings;

  while (!reader.allDone()) {
    double initTime = window[0].t;
    double endTime = initTime + cfg.coincidenceWindow;
    double timeMargin = (cfg.tRes > 0.0) ? (5.0 * cfg.tRes / constants::FWHM_TO_SIGMA) : 0.0;

    size_t oldSize = window.size();
    reader.fillWindow(window, endTime, timeMargin);

    if (window.size() > 10000) {
      size_t nSing = window.size();
      printf("UPS start!\n"
             "%s\n%s\n%s\n%s\n%s\n"
             ".\n.\n.\n"
             "%s\n%s\n%s\n%s\n%s\n"
             "UPS end!\n",
             window[0].stringify().c_str(), window[1].stringify().c_str(),
             window[2].stringify().c_str(), window[3].stringify().c_str(),
             window[4].stringify().c_str(),
             window[nSing - 5].stringify().c_str(),
             window[nSing - 4].stringify().c_str(),
             window[nSing - 3].stringify().c_str(),
             window[nSing - 2].stringify().c_str(),
             window[nSing - 1].stringify().c_str());
      return -1;
    }

    if (window.size() > oldSize) {
      std::sort(window.begin() + oldSize, window.end());
      std::inplace_merge(window.begin(), window.begin() + oldSize,
                         window.end());
    }

    double trueEndTime = window[0].t + cfg.coincidenceWindow;
    size_t windowSingles = 0;
    for (size_t i = 0; i < window.size(); ++i) {
      if (window[i].t < trueEndTime) {
        ++windowSingles;
      } else {
        break;
      }
    }

    if (windowSingles == 3) {
      ++nTriples;
    } else if (windowSingles == 4) {
      ++nQuadruples;
    } else if (windowSingles >= 5) {
      ++nQuintuplesOrMore;
    }

    if ((cfg.discardMultiples == 1 && windowSingles == 3) ||
        (cfg.discardMultiples == 2 && windowSingles >= 3)) {
      window.erase(window.begin() + 1, window.begin() + windowSingles);
      windowSingles = 1;
    }

    CoincidenceResult result = engine.findCoincidence(window, windowSingles);

    CoincidenceEvent c;
    c.time = window[0].t;
    c.amount = 1.0;
    c.gate_flag = 0;
    bool acceptedCoincidence = false;

    if (result.iCoincidence > 0) {
      // Possible CoincidenceEvent
      std::array<double, 3> p1;
      std::array<double, 3> p2;
      std::array<double, 3> lor0;
      std::array<double, 3> lor1;

      if (cfg.projectionMethod == ProjectionMethod::NONE) {
        lor0 = std::array<double, 3>{window[0].x, window[0].y, window[0].z};
        lor1 = std::array<double, 3>{window[result.iCoincidence].x,
                                     window[result.iCoincidence].y,
                                     window[result.iCoincidence].z};

        double p1Aux[3] = {window[0].x, window[0].y, window[0].z};
        p1 = toLocal(p1Aux, geo.detectorSizeX, geo.detectorSizeY,
                     cfg.detectorDepth, geo.ringRad,
                     geo.Rz[window[0].module % geo.modPerRing],
                     geo.Dz[window[0].module / geo.modPerRing], pBlur, gen);

        double p2Aux[3] = {window[result.iCoincidence].x,
                           window[result.iCoincidence].y,
                           window[result.iCoincidence].z};
        p2 =
            toLocal(p2Aux, geo.detectorSizeX, geo.detectorSizeY,
                    cfg.detectorDepth, geo.ringRad,
                    geo.Rz[window[result.iCoincidence].module % geo.modPerRing],
                    geo.Dz[window[result.iCoincidence].module / geo.modPerRing],
                    pBlur, gen);
      } else {
        double p1Blur[3] = {window[0].x + pBlur(gen), window[0].y + pBlur(gen),
                            window[0].z + pBlur(gen)};
        double p2Blur[3] = {window[result.iCoincidence].x + pBlur(gen),
                            window[result.iCoincidence].y + pBlur(gen),
                            window[result.iCoincidence].z + pBlur(gen)};

        lor0 = std::array<double, 3>{p1Blur[0], p1Blur[1], p1Blur[2]};
        lor1 = std::array<double, 3>{p2Blur[0], p2Blur[1], p2Blur[2]};

        p1 =
            project(p1Blur, p2Blur, geo.detectorSizeX, geo.detectorSizeY,
                    cfg.detectorDepth, geo.ringRad,
                    geo.Rz[window[0].module % geo.modPerRing],
                    geo.Dz[window[0].module / geo.modPerRing],
                    cfg.projectionMethod == ProjectionMethod::EXTEND_DETECTOR);
        p2 =
            project(p2Blur, p1Blur, geo.detectorSizeX, geo.detectorSizeY,
                    cfg.detectorDepth, geo.ringRad,
                    geo.Rz[window[result.iCoincidence].module % geo.modPerRing],
                    geo.Dz[window[result.iCoincidence].module / geo.modPerRing],
                    cfg.projectionMethod == ProjectionMethod::EXTEND_DETECTOR);
      }

      if (p1[0] >= geo.detectorSizeX || p1[0] <= 0.0 ||
          p1[1] >= geo.detectorSizeY || p1[1] <= 0.0 ||
          p1[2] >= cfg.detectorDepth || p1[2] < -1.0e-3 ||
          p2[0] >= geo.detectorSizeX || p2[0] <= 0.0 ||
          p2[1] >= geo.detectorSizeY || p2[1] <= 0.0 ||
          p2[2] >= cfg.detectorDepth || p2[2] < -1.0e-3) {

        if (cfg.projectionMethod != ProjectionMethod::EXTEND_DETECTOR) {
          char auxBuffer[1000];
          snprintf(
              auxBuffer, 1000,
              "Warning: Inconsistent CoincidenceEvent data.\n"
              " Original P1: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
              " Original P2: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
              "Resulting P1: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
              "Resulting P2: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
              "P1 Ring: %d, Module: %d (%d)\n"
              "P2 Ring: %d, Module: %d (%d)\n",
              window[0].e, window[0].x, window[0].y, window[0].z, window[0].t,
              window[result.iCoincidence].e, window[result.iCoincidence].x,
              window[result.iCoincidence].y, window[result.iCoincidence].z,
              window[result.iCoincidence].t, window[0].e, p1[0], p1[1], p1[2],
              window[0].t, window[result.iCoincidence].e, p2[0], p2[1], p2[2],
              window[result.iCoincidence].t, window[0].module / geo.modPerRing,
              window[0].module % geo.modPerRing, window[0].module,
              window[result.iCoincidence].module / geo.modPerRing,
              window[result.iCoincidence].module % geo.modPerRing,
              window[result.iCoincidence].module);
          warnings += auxBuffer;
        }
      } else {
        acceptedCoincidence = true;
        ++nCoincidences;

        if (window[0].hist != window[result.iCoincidence].hist) {
          ++nRandoms;
          hists.histR[result.iPair]++;
        }

        maxTimestamp = std::max(maxTimestamp, window[0].t);
        maxTimestamp = std::max(maxTimestamp, window[result.iCoincidence].t);

        unsigned normBin1X = 0, normBin1Y = 0, normBin2X = 0, normBin2Y = 0;
        size_t lorIdx = 0;

        if (pairs[result.iPair].a ==
            (cfg.useLogicalDetectors
                 ? window[0].module
                 : sim2devModule(geo.modPerRing, window[0].module))) {
          c.energy1 = window[0].e;
          c.xPosition1 = static_cast<unsigned short>(p1[0] / geo.dBinSizeX);
          c.yPosition1 = static_cast<unsigned short>(p1[1] / geo.dBinSizeY);

          c.energy2 = window[result.iCoincidence].e;
          c.xPosition2 = static_cast<unsigned short>(p2[0] / geo.dBinSizeX);
          c.yPosition2 = static_cast<unsigned short>(p2[1] / geo.dBinSizeY);

          normBin1X = static_cast<unsigned>(p1[0] / geo.dBinSizeNormX);
          normBin1Y = static_cast<unsigned>(p1[1] / geo.dBinSizeNormY);

          normBin2X = static_cast<unsigned>(p2[0] / geo.dBinSizeNormX);
          normBin2Y = static_cast<unsigned>(p2[1] / geo.dBinSizeNormY);
        } else {
          c.energy1 = window[result.iCoincidence].e;
          c.xPosition1 = static_cast<unsigned short>(p2[0] / geo.dBinSizeX);
          c.yPosition1 = static_cast<unsigned short>(p2[1] / geo.dBinSizeY);

          c.energy2 = window[0].e;
          c.xPosition2 = static_cast<unsigned short>(p1[0] / geo.dBinSizeX);
          c.yPosition2 = static_cast<unsigned short>(p1[1] / geo.dBinSizeY);

          normBin1X = static_cast<unsigned>(p2[0] / geo.dBinSizeNormX);
          normBin1Y = static_cast<unsigned>(p2[1] / geo.dBinSizeNormY);

          normBin2X = static_cast<unsigned>(p1[0] / geo.dBinSizeNormX);
          normBin2Y = static_cast<unsigned>(p1[1] / geo.dBinSizeNormY);
        }

        if (normBin1X >= cfg.nBinsNormX || normBin2X >= cfg.nBinsNormX ||
            normBin1Y >= cfg.nBinsNormY || normBin2Y >= cfg.nBinsNormY) {
          printf("Unexpected bin index\n");
          return 1;
        }

        normBin1X = cfg.nBinsNormX - 1 - normBin1X;
        normBin2X = cfg.nBinsNormX - 1 - normBin2X;

        c.pair = result.iPair;

        lorIdx = normBin1X + normBin1Y * cfg.N + normBin2X * cfg.N * cfg.N +
                 normBin2Y * cfg.N * cfg.N * cfg.N +
                 c.pair * cfg.N * cfg.N * cfg.N * cfg.N;

        hists.accumulateCoincidence(window[0], window[result.iCoincidence], p1,
                                    p2, c, lorIdx);
      }

      window.erase(window.begin() + result.iCoincidence);
      window.erase(window.begin());

    } else {
      // Single
      ++nSingles;

      double p[3] = {window[0].x, window[0].y, window[0].z};
      matmul3D(geo.Rz[window[0].module % geo.modPerRing].data(), p);
      p[2] -= geo.Dz[window[0].module / geo.modPerRing];
      p[0] -= geo.ringRad;
      p[1] += geo.detectorSizeX05;
      p[2] = geo.detectorSizeY - p[2];

      if (p[1] > geo.detectorSizeX * 1.01 || p[1] < -1.0e-3 ||
          p[2] > geo.detectorSizeY * 1.01 || p[2] < -1.0e-3 ||
          p[0] > cfg.detectorDepth * 1.01 || p[0] < -1.0e-3) {
        char auxBuffer[1000];
        snprintf(auxBuffer, 1000,
                 "Warning: Inconsistent SingleEvent data.\n"
                 " Original P: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                 "Resulting P: %15.5E %15.5E %15.5E %15.5E %15.15E\n"
                 "Ring: %d, Module: %d (%d)\n",
                 window[0].e, window[0].x, window[0].y, window[0].z,
                 window[0].t, window[0].e, p[1], p[2], p[0], window[0].t,
                 window[0].module / geo.modPerRing,
                 window[0].module % geo.modPerRing, window[0].module);
        warnings += auxBuffer;
      } else {
        c.energy1 = window[0].e;
        c.xPosition1 = static_cast<unsigned>(p[1] / geo.dBinSizeX);
        c.yPosition1 = static_cast<unsigned>(p[2] / geo.dBinSizeY);
        c.energy2 = 0.0;
        c.xPosition2 = 150;
        c.yPosition2 = 150;
        c.pair = window[0].module | 0x100;

        hists.accumulateSingle(window[0], {p[0], p[1], p[2]});
      }

      window.erase(window.begin());
    }

    if (cfg.outputFormat != OutputFormat::LM_ONLY_COINCIDENCES ||
        acceptedCoincidence) {
      writer->writeCoincidence(c);
    }

    if (nSingles + nCoincidences >= cfg.maxCounts)
      break;

    if (window.size() == 0 && !reader.allDone()) {
      window.push_back(reader.seedFirst());
    }

    if (iterations % 100000 == 0) {
      printf("\033[2J\033[1;1H");
      printf("%sSimultaneous processing singles: %zu\n%s\n",
             reader.progressString().c_str(), window.size(), warnings.c_str());
      warnings.clear();
    }
    ++iterations;
  }

  printf("%s%s\n", reader.progressString().c_str(), warnings.c_str());
  printf("Required iterations to process all files: %lu\n", iterations);
  printf("Number of coincidences: %lu\n", nCoincidences);
  printf("Number of singles     : %lu\n", nSingles);
  printf("Number of randoms     : %lu\n", nRandoms);
  printf("Number of triples     : %lu\n", nTriples);
  printf("Number of quadruples  : %lu\n", nQuadruples);
  printf("Number of 5+ multiples: %lu\n", nQuintuplesOrMore);

  writer->finalize();
  hists.writeAll("modules", cfg.header);

  return 0;
  } catch (const std::exception &e) {
    printf("\nFatal error: %s\n", e.what());
    return 1;
  }
}
