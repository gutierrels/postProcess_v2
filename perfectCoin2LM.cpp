
#include "common.hh"

int main(int argc, char** argv){

    printf("LM Header size  : %lu\n", static_cast<unsigned long>(sizeof(LMHeader)));
    printf("Coincidence size: %lu\n", static_cast<unsigned long>(sizeof(coincidence)));

    //Open and read information file
    FILE* finfo = fopen(argv[1], "r");
    if(finfo == nullptr){
        printf("Unable to open information file '%s'\n", argv[1]);
        return -2;
    }

    // + Constants

    //Detector properties
    double detectorDepth = 1.0;
    
    const double eRes = 0.06072650872059336;
    const double pRes = 0.016; //cm
    
    double coincidenceWindow = 5.0e-9;
    double energyWindow = 0.50;
    unsigned nBinsX = 300;
    unsigned nBinsY = 300;

    const bool applyNormSymmetry = false;
    const unsigned nBinsNormX = 32;
    const unsigned nBinsNormY = 32;

    const double emin = 511.0e3*(1.0-energyWindow);
    const double emax = 511.0e3*(1.0+energyWindow);

    const double eminkeV = emin/1.0e3;

    // Tallies
    unsigned long nCoincidences = 0;
    unsigned long nSingles = 0;

    //Create random gaussian distributions
    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution normDist{0.0, 1.0};
    std::normal_distribution pBlur{0.0, pRes};
    
    //LM Header
    LMHeader header;
    
    unsigned long long nCounts;
    const unsigned long long maxCounts = 10000000000000000;
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

    //Other LM Header parameters
    strcpy(header.identifier, "Simulation");
    header.rawCounts = static_cast<double>(nCounts);
    //header.rawCounts = 0.0;
    header.acqTime = 28800.0;
    header.activity = 235.7057675982394;  // In micro Curie
    strcpy(header.isotope, "FDG");
    header.startTime = 0.0;
    header.measurementTime = 0.0;
    header.isotopeHalfLife = 6586.26;
    header.weight = 1.0;
    header.maxTemp = 26.0;
    header.percentLoss = 0.0; //??
    //header.reserved = ;
    header.version[0] = 7;
    header.version[1] = 1;
    header.calibrationID = 1;
    header.gatePeriod = 0;
    header.DOILayer = 1;
    header.method = 6;
    header.StudyID = 1;

    const double detectorSizeX = header.detectorSizeX;
    const double detectorSizeY = header.detectorSizeY;
    const double detectorSizeX05 = header.detectorSizeX/2.0;
    const double detectorSizeY05 = header.detectorSizeY/2.0;
    const double dBinSizeX = header.detectorSizeX/static_cast<double>(nBinsX);
    const double dBinSizeY = header.detectorSizeY/static_cast<double>(nBinsY);
    const double dBinSizeNormX = header.detectorSizeX/static_cast<double>(nBinsNormX);
    const double dBinSizeNormY = header.detectorSizeY/static_cast<double>(nBinsNormY);    
    const double ringRad = header.detectorDistance/2.0;
    const double ringDistance = header.ringDistance;

    //Convert distances to mm
    header.detectorSizeX *= 10.0;
    header.detectorSizeY *= 10.0;
    header.detectorDistance *= 10.0;
    header.ringDistance *= 10.0;

    const unsigned nRings = header.ringNumber;
    const unsigned modPerRing = header.moduleNumber;
    const unsigned totalMods = header.moduleNumber * header.ringNumber;    

    // + Create transformations for each module and ring
    std::vector<double> Dz(nRings);
    for(unsigned iRing = 0; iRing < nRings; ++iRing){
      Dz[iRing] = iRing*(ringDistance + detectorSizeY);
    }
    
    std::vector<std::array<double,9>> Rz(modPerRing);
    constexpr double pi = 3.14159265359;
    const double angleStep = 2.0 * pi / static_cast<double>(modPerRing);    
    for(unsigned im = 0; im < modPerRing; ++im){

	//Calculate angle
	double angle =  static_cast<double>(im)*angleStep;
	double cangle = cos(angle);
	double sangle = sin(angle);

	//Create rotation matrix
	const std::array<double, 9> localRz = {
	  cangle,   -sangle, 0.0,
	  sangle,    cangle, 0.0,
	  0.0,          0.0, 1.0
	};

	Rz[im] = localRz;
    }

    header.print();

    // + Read detector pair indexes
    std::vector<detPair> pairIndexes;
    FILE* fdetPair = fopen("detectorList.txt","r");
    if(fdetPair == nullptr){
        printf("Unable to open detector pair list file 'detectorList.txt'\n");
        return -2;
    }
    
    char line[100];
    unsigned iline = 0;
    while(fgets(line,100,fdetPair) != nullptr){
      ++iline;
      unsigned a, b;
      if(sscanf(line, "%*u %u %u", &a, &b) != 2){
	printf("Corrupted detector pair file 'detectorList.txt'. "
	       "Check line %u: %s\n", iline, line);
	return -2;
      }
      pairIndexes.emplace_back(a,b);
    }
    fclose(fdetPair);
  

    // Create the LM file to be filled
    FILE* fLM = fopen("data.lm", "wb");
    if(fLM == nullptr){
      printf("Unable to create LM file 'data.lm'\n");
      return -1;
    }
    //Write header
    fwrite(&header, sizeof(LMHeader), 1, fLM);

    // Create Castor LM file to be filled
    FILE* fcLM = fopen("cdata.lm", "wb");
    if(fcLM == nullptr){
      printf("Unable to create LM file 'cdata.lm'\n");
      return -1;
    }

    //Open data file
    FILE* fin = fopen(argv[2], "rb");
    if(fin == nullptr){
      printf("Unable to open data file '%s'\n", argv[2]);
      return -2;
    }

    //Read pairs
    uint64_t time = 0;
    constexpr const size_t pairSize = 6*sizeof(float) + 2*sizeof(unsigned);
    float auxPos[6];
    unsigned auxIbody[2];
    unsigned elementsRead;
    while((elementsRead = fread(auxPos, sizeof(float), 6, fin)) > 0){

      if(elementsRead != 6){
	printf("Error: Corrupted data!\n"
	       "       expected 6 floats\n"
	       "       read: %u\n",
	       elementsRead);
	return -3;
      }

      //Read bodies
      elementsRead = fread(auxIbody, sizeof(unsigned), 2, fin);
      if(elementsRead != 2){
	printf("Error: Corrupted data!\n"
	       "       expected 2 unsigned\n"
	       "       read: %lu\n",
	       elementsRead);
	return -3;
      }
      if(auxIbody[0] >= totalMods){
	printf("Error: Body index (%u) greater than available modules (%u)\n",
	       auxIbody[0],totalMods);
	return -3;
      }
      if(auxIbody[1] >= totalMods){
	printf("Error: Body index (%u) greater than available modules (%u)\n",
	       auxIbody[1],totalMods);
	return -3;
      }

      //Project points to local coordinates
      double pos[3] = {auxPos[0], auxPos[1], auxPos[2]};
      std::array<double,3> p1 = toLocal(pos,
					detectorSizeX, detectorSizeY,
					detectorDepth, ringRad,
					Rz[auxIbody[0] % modPerRing],
					Dz[auxIbody[0] / modPerRing]);

      pos[0] = auxPos[3];
      pos[1] = auxPos[4];
      pos[2] = auxPos[5];
      std::array<double,3> p2 = toLocal(pos,
					detectorSizeX, detectorSizeY,
					detectorDepth, ringRad,
					Rz[auxIbody[1] % modPerRing],
					Dz[auxIbody[1] / modPerRing]);

      coincidenceCastor ccoin;
      ccoin.det1 = castorCrystalIndex(auxIbody[0], modPerRing,
				      nBinsNormX, nBinsNormY,
				      dBinSizeNormX, dBinSizeNormY,
				      p1);
	  
      ccoin.det2 = castorCrystalIndex(auxIbody[1], modPerRing,
				      nBinsNormX, nBinsNormY,
				      dBinSizeNormX, dBinSizeNormY,
				      p2);

      ccoin.timestamp = time;

      fwrite(&ccoin, sizeof(coincidenceCastor), 1, fcLM);

      //Write bruker coincidence
      unsigned iPair = getPair(pairIndexes,
			       sim2devModule(modPerRing,auxIbody[0]),
			       sim2devModule(modPerRing,auxIbody[1]));
      if(iPair < pairIndexes.size()){
	coincidence c;
	c.time = static_cast<double>(time);
	c.amount = 1.0;
	c.gate_flag = 0;
      
	c.energy1 = 511.0;
	c.xPosition1 = static_cast<unsigned short>(p1[0]/dBinSizeX);
	c.yPosition1 = static_cast<unsigned short>(p1[1]/dBinSizeY);
	  
	c.energy2 = 511.0;
	c.xPosition2 = static_cast<unsigned short>(p2[0]/dBinSizeX);
	c.yPosition2 = static_cast<unsigned short>(p2[1]/dBinSizeY);

	c.pair = iPair;

	//Save it
	fwrite(&c, sizeof(coincidence), 1, fLM);
      }
      
      ++time;      
    }
}
