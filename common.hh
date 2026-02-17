 
#ifndef __BRUKER_PENRED_POSTPROCESS__
#define __BRUKER_PENRED_POSTPROCESS__

#include <cstdio>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <random>
#include <thread>
#include <filesystem>
#include <iostream>
#include <cstdint>
#include <array>
#include <fstream>
#include <map>

#define WITH_HIST 1

namespace COINCIDENCE_METHOD{
  enum COINCIDENCE_METHOD : unsigned{
    TAKE_WINNER_IF_ALL_ARE_GOODS = 0,
    TAKE_CLOSEST,
    TAKE_SAME_HISTORY,
    TAKE_SAME_HISTORY_511,
  };
}

namespace OUTPUT_FORMAT{
  enum OUTPUT_FORMAT : unsigned{
    BRUKER_LM = 0,
    BRUKER_LM_ONLY_COINCIDENCES,
    CASTOR,
  };
}

namespace PROJECTION_METHOD{
  enum PROJECTION_METHOD : unsigned{
    NONE = 0,
    FIT_IN_DETECTOR,
    EXTEND_DETECTOR,
  };
}

namespace ATTENUATION_METHOD{
  enum ATTENUATION_METHOD : unsigned{
    NONE = 0,
    CYLINDER,
  };
}

namespace NORMALIZATION_METHOD{
  enum NORMALIZATION_METHOD : unsigned{
    NONE = 0,
    CREATE,
    APPLY,
  };
}

inline unsigned sim2devModule(const unsigned ringModules, const unsigned imod){
  const unsigned iRing = imod/ringModules;
  return iRing*ringModules + ((imod + (ringModules-2)) % ringModules);
}

inline unsigned castorCrystalIndex(const unsigned iModuleGlob,
				   const unsigned modulesPerRing,
				   const double nxBins, const double nyBins,
				   const double dxBin, const double dyBin,
				   const std::array<double, 3>& pLocal){
  const unsigned cristalsPerRing = nxBins*modulesPerRing;
  const unsigned iModuleDev = sim2devModule(modulesPerRing, iModuleGlob);
  //const unsigned iModuleDev = iModuleGlob;
  const unsigned iRing = iModuleDev / modulesPerRing;
  const unsigned iModule = iModuleDev % modulesPerRing;

  unsigned ix = nxBins-1 - static_cast<unsigned>(pLocal[0]/dxBin);
  unsigned iy = nyBins-1 - static_cast<unsigned>(pLocal[1]/dyBin);

  unsigned iCastorRing = iRing*nyBins + iy;
  return iCastorRing*cristalsPerRing + iModule*nxBins + ix;
}

inline std::pair<uint32_t, uint32_t> castorLOR(const unsigned iModule1Glob,
					       const unsigned iModule2Glob,
					       const std::array<double, 3>& p1Local,
					       const std::array<double, 3>& p2Local,
					       const unsigned modulesPerRing,
					       const double nxBins, const double nyBins,
					       const double dxBin, const double dyBin){
  
  unsigned det1 = castorCrystalIndex(iModule1Glob,modulesPerRing,
				     nxBins,nyBins,dxBin,dyBin,p1Local);
  unsigned det2 = castorCrystalIndex(iModule2Glob,modulesPerRing,
				     nxBins,nyBins,dxBin,dyBin,p2Local);

  return std::pair<uint32_t,uint32_t>(det1, det2);
}

struct single{

  static constexpr const size_t bufferSize = 5*sizeof(float) + sizeof(double) + 3*sizeof(uint8_t) + sizeof(unsigned long long);
  
  float e, x, y, z;
  double t;
  std::array<uint8_t, 3> info;
  unsigned long long hist;
  unsigned module;

  single() : t(1.0e35){}

  inline int readSingle(FILE* f, const unsigned m, const double emin, const double emax){
    //Reset time
    t = 1.0e35;
    //Reset energy
    e = -1.0;
    //Set module
    module = m;

    while(e < emin || e > emax){
      //Read next single

      fread(&e, sizeof(float), 1, f);
      fread(&x, sizeof(float), 1, f);
      fread(&y, sizeof(float), 1, f);
      fread(&z, sizeof(float), 1, f);

      //Skip weight
      float w;
      fread(&w, sizeof(float), 1, f);

      fread(&t, sizeof(double), 1, f);
      fread(info.data(), sizeof(uint8_t), 3, f);
      
      if(fread(&hist, sizeof(unsigned long long), 1, f) != 1){
	printf("End of file reached\n");
	fflush(stdout);
	t = 1.0e35;
	return 1;	
      }
    }
    return 0;
  }

  inline int readPenRed(FILE* f, const unsigned m, const double emin, const double emax,
			const double eRes, std::normal_distribution<double>& normDist, std::mt19937& gen){
    //Reset time
    t = 1.0e35;
    //Reset energy
    e = -1.0;
    //Set module
    module = m;

    char buffer[bufferSize];
    
    while(e < emin || e > emax){

      if(fread(buffer, 1, bufferSize, f) == bufferSize){

	size_t pos = 0;
		
	memcpy(&e, buffer, sizeof(float));
	pos += sizeof(float);
	memcpy(&x, buffer+pos, sizeof(float));
	pos += sizeof(float);

	memcpy(&y, buffer+pos, sizeof(float));
	pos += sizeof(float);
		
	memcpy(&z, buffer+pos, sizeof(float));
	pos += sizeof(float);

	//Skip weight
	//memcpy(&w, buffer+pos, sizeof(float)); 
	//pos += sizeof(float);
		
	memcpy(&t, buffer+pos, sizeof(double));
	pos += sizeof(double);

	//memcpy(info.data(), buffer+pos, 3*sizeof(uint8_t));
	//pos += 3*sizeof(uint8_t);
	
	memcpy(&hist, buffer+pos, sizeof(unsigned long long));
      }
      else{
	//printf("End of file reached\n");
	//fflush(stdout);
	t = 1.0e35;
	return 1;	
      }

      //Apply energy bluring
      double sigma = e * eRes / 2.355;
      double normVal = normDist(gen);
      e += sigma*normVal;
    }

    //Convert energy to keV
    e /= 1.0e3;
    
    return 0;
  }

  inline std::string stringify() const {

    char data[300];
    snprintf(data, 300, "%15.5E %15.5E %15.5E "
	     "%15.5E %25.15E %u %u %u %llu\n",
	     e, x, y, z, t,
	     info[0], info[1], info[2], hist);
    return std::string(data);
  }
  
  inline bool operator>(const single& o) const {
    return t > o.t;
  }
  inline bool operator<(const single& o) const {
    return t < o.t;
  }
  inline bool operator==(const single& o) const {
    return t == o.t;
  }
  
};

#pragma pack(push)
#pragma pack(4)
struct coincidenceCastor{
  uint32_t timestamp; //ms
  uint32_t det1, det2;
};
#pragma pack(pop)

struct coincidence{
  double         time;
  float          energy1;
  float          energy2;
  float          amount;
  unsigned short xPosition1;
  unsigned short yPosition1;
  unsigned short xPosition2;
  unsigned short yPosition2;
  unsigned short pair;
  unsigned short gate_flag;

  inline void write(FILE* fout) const {
    fprintf(fout, "%03u %.3E %.15E "
	    "%03u %03u %.5E "
	    "%03u %03u %.5E "
	    "%03u\n",
	    pair, amount, time,
	    xPosition1, yPosition1, energy1,
	    xPosition2, yPosition2, energy2,
	    gate_flag);
  }

  inline void print() const{
    printf("%03u %.3E %.15E "
	   "%03u %03u %.5E "
	   "%03u %03u %.5E "
	   "%03u\n",
	   pair, amount, time,
	   xPosition1, yPosition1, energy1,
	   xPosition2, yPosition2, energy2,
	   gate_flag);    
  }

  static inline bool read(coincidence& c, FILE* fin){
    return fread(&c, sizeof(coincidence), 1, fin) == sizeof(coincidence);
  }
};

#pragma pack(push)
#pragma pack(4)
struct LMHeader
{
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
  char unused_1;  //padding
  double gatePeriod;
  short DOILayer;
  short method;
  short StudyID;
  char unused_2[6];


  static inline bool read(LMHeader& header, FILE* fin){
    return fread(&header, sizeof(LMHeader), 1, fin) == sizeof(LMHeader);
  }

  void print() const{
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
	   identifier, rawCounts, acqTime, activity, isotope,
	   detectorSizeX, detectorSizeY, startTime, measurementTime,
	   moduleNumber, ringNumber, ringDistance, detectorDistance,
	   isotopeHalfLife, weight, maxTemp, percentLoss, reserved,
	   version, calibrationID, gatePeriod, DOILayer, method, StudyID);
  }
};
#pragma pack(pop)

struct detPair {
  unsigned a,b;

  detPair() : a(0), b(0){}
  detPair(const unsigned _a, const unsigned _b) : a(_a), b(_b){}

  inline unsigned other(const unsigned c) const {
    if(a == c)
      return b;
    if(b == c)
      return a;
    return c;
  }
  inline bool operator==(const detPair& o) const {
    return (a == o.a && b == o.b) || (a == o.b && b == o.a);
  }
};

inline unsigned getPair(const std::vector<detPair>& list, const unsigned a, const unsigned b){

  const detPair p(a,b);
  for(unsigned i = 0; i < list.size(); ++i){
    if(p == list[i])
      return i;
  }
  return list.size();
}

std::array<double,3> project(const double* p1Orig, const double* p2Orig,
			     const double detSizeX, const double detSizeY,
			     const double detDepth,
			     const double ringRad,
			     const std::array<double, 9>& Rz, double Dz,
			     const bool extendDetector);

std::array<double,3> toLocal(const double* p1Orig,
			     const double detSizeX, const double detSizeY,
			     const double detDepth,
			     const double ringRad,
			     const std::array<double, 9>& Rz, double Dz);
inline std::array<double,3> toLocal(const double* p1Orig,
				    const double detSizeX, const double detSizeY,
				    const double detDepth,
				    const double ringRad,
				    const std::array<double, 9>& Rz, double Dz,
				    std::normal_distribution<double>& pBlur, std::mt19937& gen){

  const std::array<double,3> local = toLocal(p1Orig, detSizeX, detSizeY, detDepth, ringRad, Rz, Dz);
  
  // Apply bluring
  unsigned count = 0;
  std::array<double,3> res;
  do{
    res[0] = local[0] + pBlur(gen);
  } while((res[0] <= 0.0 || res[0] >= detSizeX) && count++ < 1000);
  do{
    res[1] = local[1] + pBlur(gen);
  } while((res[1] <= 0.0 || res[1] >= detSizeY) && count++ < 1000);
  do{
    res[2] = local[2] + pBlur(gen);
  } while((res[2] <= 0.0 || res[2] >= detDepth) && count++ < 1000);

  if(count >= 1000)
    printf("Error bluring!\n");
  
  return res;
}

inline std::string stringifyFileProgress(FILE* f, const std::uintmax_t fSize,
					 const char* prefix){

  std::stringstream ss;

  ss << std::fixed << std::setprecision(2);
  
  ss << prefix << " [";
	
  //Calculate file progress
  double progress = 1.0;
  if(f != nullptr){
    progress = static_cast<double>(ftell(f))/static_cast<double>(fSize);
  }

  constexpr const int ibarWidth = 50;
  constexpr const double barWidth = static_cast<double>(ibarWidth);
  
  int pos = static_cast<int>(barWidth*progress);
  for(int ib = 0; ib < pos; ++ib){
    ss << "=";
  }
  if(pos < ibarWidth-1){
    ss << ">";
    for(int ib = pos; ib < ibarWidth; ++ib){
      ss << " ";
    }
  }
  ss << "] " << progress*100.0 << " %\n";

  return ss.str();
}

inline void matmul3D(const double A[9], double B[3]){

  double C[3];

  for(int j = 0; j < 3; j++){
    C[j] = 0.0;
    for(int i = 0; i < 3; i++){
      C[j] += A[3*j+i]*B[i];
    }
  }
  for(int i = 0; i < 3; i++)
    B[i] = C[i];

}


inline double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline std::array<double, 3> subtract(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

inline std::array<double, 3> add(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

inline std::array<double, 3> scale(const std::array<double, 3>& v, double s) {
    return {v[0]*s, v[1]*s, v[2]*s};
}

double attenuationFactorWithAir(const std::array<double, 3>& p1,
				const std::array<double, 3>& p2,
				double mu_cylinder,
				double mu_air,
				const std::array<double, 3>& cylinder_origin,
				double radius,
				double height);

#endif
