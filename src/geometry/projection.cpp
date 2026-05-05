
#include "geometry/projection.hh"
#include "common/math_utils.hh"
#include <algorithm>
#include <cmath>

std::array<double, 3> project(const double *p1Orig, const double *p2Orig,
                              const double detSizeX, const double detSizeY,
                              const double detDepth, const double ringRad,
                              const std::array<double, 9> &Rz, const double Dz,
                              const bool extendDetector) {

  const double detSizeX05 = detSizeX / 2.0;

  double p1[3] = {p1Orig[0], p1Orig[1], p1Orig[2]};
  double p2[3] = {p2Orig[0], p2Orig[1], p2Orig[2]};

  // Rotate both positions to locate the first point within the first ring
  // detector
  matmul3D(Rz.data(), p1);
  matmul3D(Rz.data(), p2);

  // Move both points to locate the first one on the first ring
  p1[2] -= Dz;
  p2[2] -= Dz;

  // Clamp position to fit detector size (Can fall outside due a previous
  // applied bluring)
  p1[0] =
      std::clamp(p1[0], ringRad * 1.000001, (ringRad + detDepth) * 0.999999);
  p1[1] = std::clamp(p1[1], -detSizeX05 * 0.999999, detSizeX05 * 0.999999);
  p1[2] = std::clamp(p1[2], 1.0e-6, detSizeY * 0.999999);

  // Calculate point distances
  std::array<double, 3> dp = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};

  double tmin = 1.0e35;
  if (extendDetector) {
    // the projection can extend beyond the detector borders
    if (dp[0] != 0.0)
      tmin = (ringRad - p1[0]) / dp[0]; // -X
  } else {
    // The projection have to fit inside the detector borders

    // Calculate the points where the line cuts the detector in each direction
    //(x=ringRad or y=+-halfDetSizeX or z=+-halfDetSizeY)
    //  x = x0 + (xf-x0)*t
    //  y = y0 + (yf-y0)*t
    //  z = z0 + (zf-z0)*t
    const std::array<double, 5> ts = {
        dp[0] == 0 ? 1.0e35
                   : (ringRad - p1[0]) / dp[0], // -X (+X is not required)
        dp[1] == 0 ? 1.0e35 : (-detSizeX05 - p1[1]) / dp[1], // -Y
        dp[1] == 0 ? 1.0e35 : (detSizeX05 - p1[1]) / dp[1],  // +Y
        dp[2] == 0 ? 1.0e35 : (-p1[2]) / dp[2],              // -Z
        dp[2] == 0 ? 1.0e35 : (detSizeY - p1[2]) / dp[2],    // +Z
    };

    for (double taux : ts) {
      if (taux >= 0.0) {
        tmin = std::min(tmin, taux);
      } else if (taux > -1.0e-6) {
        tmin = 0.0;
        break;
      }
    }
  }

  const double t = tmin < 1.0e30 ? tmin : 0.0;

  // Project the point
  p1[0] = p1[0] + dp[0] * t * 0.999; // X
  p1[1] = p1[1] + dp[1] * t * 0.999; // Y
  p1[2] = p1[2] + dp[2] * t * 0.999; // Z

  // Transform the point to local detector coordinates
  p1[0] -= ringRad;
  p1[1] += detSizeX05;
  p1[2] = detSizeY - p1[2];

  // Return final point in local coordinates
  // Note: In local coordinates (x,y,DOI) = (p[1],p[2],p[0])
  return {p1[1], p1[2], p1[0]};
}

std::array<double, 3> toLocal(const double *p1Orig, const double detSizeX,
                              const double detSizeY, const double detDepth,
                              const double ringRad,
                              const std::array<double, 9> &Rz, double Dz) {

  const double detSizeX05 = detSizeX / 2.0;

  double p1[3] = {p1Orig[0], p1Orig[1], p1Orig[2]};

  // Rotate both positions to locate the first point within the first ring
  // detector
  matmul3D(Rz.data(), p1);

  // Move both points to locate the first one on the first ring
  p1[2] -= Dz;

  // Check bounds
  if (p1[0] > (ringRad + detDepth) * 1.001 || p1[0] < ringRad * 0.999 ||
      p1[1] > detSizeX05 * 1.001 || p1[1] < -detSizeX05 * 1.001 ||
      p1[2] > detSizeY * 1.001 || p1[2] < -1.0e-3) {
    printf("Error: Unconsistent coincidence data.\n");
    printf(" Original P1: %15.5E %15.5E %15.5E\n", p1Orig[0], p1Orig[1],
           p1Orig[2]);
    printf("Resulting P1: %15.5E %15.5E %15.5E\n", p1[0], p1[1], p1[2]);
    printf("Valid transformed P1 intervals: ( [%15.5E,%15.5E], "
           "[%15.5E,%15.5E], [%15.5E,%15.5E] )\n",
           ringRad, ringRad + detDepth, -detSizeX05, detSizeX05, 0.0, detSizeY);
  }

  // Clamp position
  p1[0] =
      std::clamp(p1[0], ringRad * 1.000001, (ringRad + detDepth) * 0.999999);
  p1[1] = std::clamp(p1[1], -detSizeX05 * 0.999999, detSizeX05 * 0.999999);
  p1[2] = std::clamp(p1[2], 1.0e-6, detSizeY * 0.999999);

  // Transform the point to local detector coordinates
  p1[0] -= ringRad;
  p1[1] += detSizeX05;
  p1[2] = detSizeY - p1[2];

  // Return final point in local coordinates
  // Note: In local coordinates (x,y,DOI) = (p[1],p[2],p[0])
  return {p1[1], p1[2], p1[0]};
}

std::array<double, 3>
toLocal(const double *p1Orig, const double detSizeX, const double detSizeY,
        const double detDepth, const double ringRad,
        const std::array<double, 9> &Rz, double Dz,
        std::normal_distribution<double> &pBlur, std::mt19937 &gen) {

  const std::array<double, 3> local =
      toLocal(p1Orig, detSizeX, detSizeY, detDepth, ringRad, Rz, Dz);

  // Apply bluring
  unsigned count = 0;
  std::array<double, 3> res;
  do {
    res[0] = local[0] + pBlur(gen);
  } while ((res[0] <= 0.0 || res[0] >= detSizeX) && count++ < 1000);
  do {
    res[1] = local[1] + pBlur(gen);
  } while ((res[1] <= 0.0 || res[1] >= detSizeY) && count++ < 1000);
  do {
    res[2] = local[2] + pBlur(gen);
  } while ((res[2] <= 0.0 || res[2] >= detDepth) && count++ < 1000);

  if (count >= 1000)
    printf("Error bluring!\n");

  return res;
}



void toLogical(single &s, const unsigned modPerRing,
               const double logicalDetSizeY,
               const std::vector<std::array<double, 9>> &Rz,
               const std::vector<double> &Dz) {

  unsigned physRing = s.module / modPerRing;
  unsigned modInRing = s.module % modPerRing;

  double p[3] = {s.x, s.y, s.z};
  matmul3D(Rz[modInRing].data(), p); // Rotar al módulo 0
  p[2] -= Dz[physRing];              // Quitar traslación axial del anillo

  bool isUpper = (p[2] >= logicalDetSizeY);

  unsigned logRing = physRing * 2 + (isUpper ? 1 : 0);
  s.module = logRing * modPerRing + modInRing;
}
