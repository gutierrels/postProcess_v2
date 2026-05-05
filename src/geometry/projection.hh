
#ifndef POSTPROCESS_GEOMETRY_PROJECTION_HH
#define POSTPROCESS_GEOMETRY_PROJECTION_HH

#include "common/types.hh"
#include <array>
#include <vector>
#include <random>

std::array<double, 3> project(const double *p1Orig, const double *p2Orig,
                              const double detSizeX, const double detSizeY,
                              const double detDepth, const double ringRad,
                              const std::array<double, 9> &Rz, double Dz,
                              const bool extendDetector);

std::array<double, 3> toLocal(const double *p1Orig, const double detSizeX,
                              const double detSizeY, const double detDepth,
                              const double ringRad,
                              const std::array<double, 9> &Rz, double Dz);

std::array<double, 3>
toLocal(const double *p1Orig, const double detSizeX, const double detSizeY,
        const double detDepth, const double ringRad,
        const std::array<double, 9> &Rz, double Dz,
        std::normal_distribution<double> &pBlur, std::mt19937 &gen);



void toLogical(single &s, const unsigned modPerRing,
               const double logicalDetSizeY,
               const std::vector<std::array<double, 9>> &Rz,
               const std::vector<double> &Dz);

#endif
