
#include "io/pair_list.hh"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

std::vector<DetectorPair> readPairList(const std::string &filename) {
  std::vector<DetectorPair> pairIndexes;
  FILE *fdetPair = fopen(filename.c_str(), "r");
  if (fdetPair == nullptr) {
    throw std::runtime_error("Unable to open detector pair list file '" + filename + "'");
  }

  printf("Reading detector pair list file '%s'\n", filename.c_str());

  char line[100];
  unsigned iline = 0;
  while (fgets(line, 100, fdetPair) != nullptr) {
    ++iline;
    unsigned a, b;
    if (sscanf(line, "%*u %u %u", &a, &b) != 2) {
      throw std::runtime_error("Corrupted detector pair file '" + filename + "'. Check line " + std::to_string(iline) + ": " + std::string(line));
    }
    pairIndexes.emplace_back(a, b);
  }
  fclose(fdetPair);
  printf("Loaded %lu detector pairs\n", pairIndexes.size());

  return pairIndexes;
}
