
#ifndef POSTPROCESS_COMMON_TYPES_HH
#define POSTPROCESS_COMMON_TYPES_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

namespace COINCIDENCE_METHOD {
enum COINCIDENCE_METHOD : unsigned {
  TAKE_WINNER_IF_ALL_ARE_GOODS = 0,
  TAKE_CLOSEST,
  TAKE_SAME_HISTORY,
  TAKE_SAME_HISTORY_511,
};
}

namespace OUTPUT_FORMAT {
enum OUTPUT_FORMAT : unsigned {
  BRUKER_LM = 0,
  BRUKER_LM_ONLY_COINCIDENCES,
};
}

namespace PROJECTION_METHOD {
enum PROJECTION_METHOD : unsigned {
  NONE = 0,
  FIT_IN_DETECTOR,
  EXTEND_DETECTOR,
};
}

namespace GENERATE_HISTOGRAM {
enum GENERATE_HISTOGRAM : unsigned {
  NONE = 0,
  LOR_INDEX,
};
}

// ---------------------------------------------------------------------------
// Module-index helpers
// ---------------------------------------------------------------------------

inline unsigned sim2devModule(const unsigned ringModules, const unsigned imod) {
  const unsigned iRing = imod / ringModules;
  return iRing * ringModules + ((imod + (ringModules - 2)) % ringModules);
}

// ---------------------------------------------------------------------------
// single
// ---------------------------------------------------------------------------

struct single {

  // PenRed Save Flags Configuration
  // Adjust these to match tallies/singles/save/* in your PenRed config
  static constexpr bool saveWeight = true;
  static constexpr bool saveMetadata = true;
  static constexpr bool saveHistory = true;

  static constexpr const size_t bufferSize =
      5 * sizeof(float) + sizeof(double) + 3 * sizeof(uint8_t) +
      sizeof(unsigned long long);

  float e, x, y, z;
  double t;
  double true_t;
  std::array<uint8_t, 3> info;
  unsigned long long hist;
  unsigned module;

  single() : t(1.0e35), true_t(1.0e35) {}

  inline int readPenRed(FILE *f, const unsigned m, const double emin,
                        const double emax, const double eRes, const double tRes,
                        std::normal_distribution<double> &normDist,
                        std::mt19937 &gen) {
    // Reset time
    t = 1.0e35;
    true_t = 1.0e35;
    // Reset energy
    e = -1.0;
    // Set module
    module = m;

    char buffer[bufferSize];

    while (e < emin || e > emax) {

      if (fread(buffer, 1, bufferSize, f) == bufferSize) {

        size_t pos = 0;

        memcpy(&e, buffer + pos, sizeof(float));
        pos += sizeof(float);
        memcpy(&x, buffer + pos, sizeof(float));
        pos += sizeof(float);
        memcpy(&y, buffer + pos, sizeof(float));
        pos += sizeof(float);
        memcpy(&z, buffer + pos, sizeof(float));
        pos += sizeof(float);

        if constexpr (saveWeight) {
          // If you ever need to store weight in the struct, do it here
          // memcpy(&w, buffer + pos, sizeof(float));
          pos += sizeof(float);
        }

        memcpy(&t, buffer + pos, sizeof(double));
        pos += sizeof(double);

        if constexpr (saveMetadata) {
          memcpy(info.data(), buffer + pos, 3 * sizeof(uint8_t));
          pos += 3 * sizeof(uint8_t);
        }

        if constexpr (saveHistory) {
          memcpy(&hist, buffer + pos, sizeof(unsigned long long));
          pos += sizeof(unsigned long long);
        }
      } else {
        // printf("End of file reached\n");
        // fflush(stdout);
        t = 1.0e35;
        true_t = 1.0e35;
        return 1;
      }

      true_t = t;

      // Apply energy bluring
      double sigma = e * eRes / 2.355;
      double normVal = normDist(gen);
      e += sigma * normVal;

      // Apply time bluring
      if (tRes > 0.0) {
        double sigma_t = tRes / 2.355;
        t += sigma_t * normDist(gen);
        if (t < 0.0)
          t = 0.0;
      }
    }

    // Convert energy to keV
    e /= 1.0e3;

    return 0;
  }

  inline std::string stringify() const {

    char data[300];
    snprintf(data, 300,
             "%15.5E %15.5E %15.5E "
             "%15.5E %25.15E %u %u %u %llu\n",
             e, x, y, z, t, info[0], info[1], info[2], hist);
    return std::string(data);
  }

  inline bool operator>(const single &o) const { return t > o.t; }
  inline bool operator<(const single &o) const { return t < o.t; }
  inline bool operator==(const single &o) const { return t == o.t; }
};

// ---------------------------------------------------------------------------
// coincidence
// ---------------------------------------------------------------------------

struct coincidence {
  double time;
  float energy1;
  float energy2;
  float amount;
  unsigned short xPosition1;
  unsigned short yPosition1;
  unsigned short xPosition2;
  unsigned short yPosition2;
  unsigned short pair;
  unsigned short gate_flag;

  inline void write(FILE *fout) const {
    fprintf(fout,
            "%03u %.3E %.15E "
            "%03u %03u %.5E "
            "%03u %03u %.5E "
            "%03u\n",
            pair, amount, time, xPosition1, yPosition1, energy1, xPosition2,
            yPosition2, energy2, gate_flag);
  }

  inline void print() const {
    printf("%03u %.3E %.15E "
           "%03u %03u %.5E "
           "%03u %03u %.5E "
           "%03u\n",
           pair, amount, time, xPosition1, yPosition1, energy1, xPosition2,
           yPosition2, energy2, gate_flag);
  }

  static inline bool read(coincidence &c, FILE *fin) {
    return fread(&c, sizeof(coincidence), 1, fin) == sizeof(coincidence);
  }
};

// ---------------------------------------------------------------------------
// LMHeader
// ---------------------------------------------------------------------------

#pragma pack(push)
#pragma pack(4)
struct LMHeader {
  char identifier[16];
  double rawCounts;
  double acqTime;
  double activity;
  char isotope[16];
  double detectorSizeX;
  double detectorSizeY;
  double startTime;
  double measurementTime;
  int moduleNumber;
  int ringNumber;
  double ringDistance;
  double detectorDistance;
  double isotopeHalfLife;
  float weight;
  float maxTemp;
  float percentLoss;
  float reserved[5];
  unsigned char version[2];
  unsigned char calibrationID;
  char unused_1; // padding
  double gatePeriod;
  short DOILayer;
  short method;
  short StudyID;
  char unused_2[6];

  static inline bool read(LMHeader &header, FILE *fin) {
    return fread(&header, sizeof(LMHeader), 1, fin) == sizeof(LMHeader);
  }

  void print() const {
    printf("Identifier        : %s\n"
           "Raw Counts        : %E\n"
           "Acquisition Time  : %6f s\n"
           "Activity          : %6f  muC\n"
           "Isotope           : %s\n"
           "Detector Size (x) : %6f mm\n"
           "Detector Size (y) : %6f mm\n"
           "Start time        : %6f s\n"
           "Measurement time  : %6f s\n"
           "Modules           : %4d\n"
           "Rings             : %4d\n"
           "Ring Distance     : %6f mm\n"
           "Detector Distance : %6f mm\n"
           "Isotope HalfLife  : %6f s\n"
           "Weight            : %6f\n"
           "Max Temperature   : %6f ºC\n"
           "Loss              : %6f %%\n"
           "Version           : %s\n"
           "Calibration ID    : %c\n"
           "Gate Period       : %6f s\n"
           "DOI Layer         : %4d\n"
           "Method            : %4d\n"
           "Study ID          : %4d\n",
           identifier, rawCounts, acqTime, activity, isotope, detectorSizeX,
           detectorSizeY, startTime, measurementTime, moduleNumber, ringNumber,
           ringDistance, detectorDistance, isotopeHalfLife, weight, maxTemp,
           percentLoss, version, calibrationID, gatePeriod, DOILayer, method,
           StudyID);
  }
};
#pragma pack(pop)

// ---------------------------------------------------------------------------
// detPair
// ---------------------------------------------------------------------------

struct detPair {
  unsigned a, b;

  detPair() : a(0), b(0) {}
  detPair(const unsigned _a, const unsigned _b) : a(_a), b(_b) {}

  inline unsigned other(const unsigned c) const {
    if (a == c)
      return b;
    if (b == c)
      return a;
    return c;
  }
  inline bool operator==(const detPair &o) const {
    return (a == o.a && b == o.b) || (a == o.b && b == o.a);
  }
};

inline unsigned getPair(const std::vector<detPair> &list, const unsigned a,
                        const unsigned b) {

  const detPair p(a, b);
  for (unsigned i = 0; i < list.size(); ++i) {
    if (p == list[i])
      return i;
  }
  return list.size();
}

#endif
