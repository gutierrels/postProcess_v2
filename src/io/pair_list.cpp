
#include "io/pair_list.hh"
#include <cstdio>
#include <cstdlib>

std::vector<detPair> readPairList(const std::string &filename) {
  std::vector<detPair> pairIndexes;
  FILE *fdetPair = fopen(filename.c_str(), "r");
  if (fdetPair == nullptr) {
    printf("Unable to open detector pair list file '%s'\n", filename.c_str());
    std::exit(-2);
  }

  printf("Reading detector pair list file '%s'\n", filename.c_str());

  char line[100];
  unsigned iline = 0;
  while (fgets(line, 100, fdetPair) != nullptr) {
    ++iline;
    unsigned a, b;
    if (sscanf(line, "%*u %u %u", &a, &b) != 2) {
      printf("Corrupted detector pair file '%s'. "
             "Check line %u: %s\n",
             filename.c_str(), iline, line);
      std::exit(-2);
    }
    pairIndexes.emplace_back(a, b);
  }
  fclose(fdetPair);
  printf("Loaded %lu detector pairs\n", pairIndexes.size());

  return pairIndexes;
}
