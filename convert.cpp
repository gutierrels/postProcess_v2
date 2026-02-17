 

#include "common.hh"

int main(int argc, char** argv){

  //Open and read information file
  FILE* finfo = fopen(argv[1], "r");
  if(finfo == nullptr){
    printf("Unable to open information file '%s'\n", argv[1]);
    return -2;
  }

  unsigned long long nCounts;
  if(fscanf(finfo,
	    "%llu %lE %lE "
	    "%lE %d %d "
	    "%lE",
	    &nCounts, &header.detectorSizeX, &header.detectorSizeY,
	    &header.detectorDistance, &header.moduleNumber, &header.ringNumber,
	    &header.ringDistance) != 7){
    printf("Corrupted information file '%s'\n", argv[1]);
    return -2;
  }
  
  
  return 0;
}
