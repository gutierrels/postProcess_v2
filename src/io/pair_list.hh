
#ifndef POSTPROCESS_IO_PAIR_LIST_HH
#define POSTPROCESS_IO_PAIR_LIST_HH

#include "common/types.hh"
#include <string>
#include <vector>

std::vector<DetectorPair> readPairList(const std::string &filename);

#endif
