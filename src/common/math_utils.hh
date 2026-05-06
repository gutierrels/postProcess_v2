
#ifndef POSTPROCESS_COMMON_MATH_UTILS_HH
#define POSTPROCESS_COMMON_MATH_UTILS_HH

#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>

inline void matmul3D(const double A[9], double B[3]) {

  double C[3];

  for (int j = 0; j < 3; j++) {
    C[j] = 0.0;
    for (int i = 0; i < 3; i++) {
      C[j] += A[3 * j + i] * B[i];
    }
  }
  for (int i = 0; i < 3; i++)
    B[i] = C[i];
}

inline double dot(const std::array<double, 3> &a,
                  const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline std::array<double, 3> subtract(const std::array<double, 3> &a,
                                      const std::array<double, 3> &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline std::array<double, 3> add(const std::array<double, 3> &a,
                                 const std::array<double, 3> &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline std::array<double, 3> scale(const std::array<double, 3> &v, double s) {
  return {v[0] * s, v[1] * s, v[2] * s};
}

inline std::string stringifyFileProgress(FILE *f, const std::uintmax_t fSize,
                                         const char *prefix) {

  std::stringstream ss;

  ss << std::fixed << std::setprecision(2);

  ss << prefix << " [";

  // Calculate file progress
  double progress = 1.0;
  if (f != nullptr) {
    progress = static_cast<double>(ftell(f)) / static_cast<double>(fSize);
  }

  constexpr const int ibarWidth = 50;
  constexpr const double barWidth = static_cast<double>(ibarWidth);

  int pos = static_cast<int>(barWidth * progress);
  for (int ib = 0; ib < pos; ++ib) {
    ss << "=";
  }
  if (pos < ibarWidth - 1) {
    ss << ">";
    for (int ib = pos; ib < ibarWidth; ++ib) {
      ss << " ";
    }
  }
  ss << "] " << progress * 100.0 << " %\n";

  return ss.str();
}

#endif
