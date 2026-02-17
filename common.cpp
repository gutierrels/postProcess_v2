#include "common.hh"

std::array<double,3> project(const double* p1Orig, const double* p2Orig,
			     const double detSizeX, const double detSizeY,
			     const double detDepth,
			     const double ringRad,
			     const std::array<double, 9>& Rz, const double Dz,
			     const bool extendDetector){

  const double detSizeX05 = detSizeX/2.0;

  double p1[3] = {p1Orig[0], p1Orig[1], p1Orig[2]};
  double p2[3] = {p2Orig[0], p2Orig[1], p2Orig[2]};
  
  //Rotate both positions to locate the first point within the first ring detector
  matmul3D(Rz.data(), p1);
  matmul3D(Rz.data(), p2);

  //Move both points to locate the first one on the first ring
  p1[2] -= Dz;
  p2[2] -= Dz;
  
  //Clamp position to fit detector size (Can fall outside due a previous applied bluring)
  p1[0] = std::clamp(p1[0],     ringRad*1.000001, (ringRad+detDepth)*0.999999);
  p1[1] = std::clamp(p1[1], -detSizeX05*0.999999,         detSizeX05*0.999999);
  p1[2] = std::clamp(p1[2],               1.0e-6,           detSizeY*0.999999);

  //Calculate point distances
  std::array<double,3> dp = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};

  double tmin = 1.0e35;
  if(extendDetector){
    //the projection can extend beyond the detector borders
    if(dp[0] != 0.0)
      tmin = (ringRad - p1[0]) / dp[0]; // -X
  }
  else{
    //The projection have to fit inside the detector borders
  
    //Calculate the points where the line cuts the detector in each direction
    //(x=ringRad or y=+-halfDetSizeX or z=+-halfDetSizeY)
    // x = x0 + (xf-x0)*t
    // y = y0 + (yf-y0)*t
    // z = z0 + (zf-z0)*t
    const std::array<double, 5> ts = {
      dp[0] == 0 ? 1.0e35 : (    ringRad - p1[0]) / dp[0], // -X (+X is not required)
      dp[1] == 0 ? 1.0e35 : (-detSizeX05 - p1[1]) / dp[1], // -Y
      dp[1] == 0 ? 1.0e35 : ( detSizeX05 - p1[1]) / dp[1], // +Y
      dp[2] == 0 ? 1.0e35 : (            - p1[2]) / dp[2], // -Z
      dp[2] == 0 ? 1.0e35 : (   detSizeY - p1[2]) / dp[2], // +Z
    };
  
    for(double taux : ts){
      if(taux >= 0.0){
	tmin = std::min(tmin,taux);
      }
      else if(taux > -1.0e-6){
	tmin = 0.0;
	break;
      }
    }
  }
  
  const double t = tmin < 1.0e30 ? tmin : 0.0;
  
  //Project the point
  p1[0] = p1[0] + dp[0]*t*0.999; //X
  p1[1] = p1[1] + dp[1]*t*0.999; //Y
  p1[2] = p1[2] + dp[2]*t*0.999; //Z

  //Transform the point to local detector coordinates
  p1[0] -= ringRad;
  p1[1] += detSizeX05;
  p1[2] = detSizeY - p1[2];
  
  //Return final point in local coordinates
  //Note: In local coordinates (x,y,DOI) = (p[1],p[2],p[0])
  return {p1[1], p1[2], p1[0]};
}

std::array<double,3> toLocal(const double* p1Orig,
			     const double detSizeX, const double detSizeY,
			     const double detDepth,
			     const double ringRad,
			     const std::array<double, 9>& Rz, double Dz){

  const double detSizeX05 = detSizeX/2.0;

  double p1[3] = {p1Orig[0], p1Orig[1], p1Orig[2]};

  //Rotate both positions to locate the first point within the first ring detector
  matmul3D(Rz.data(), p1);

  //Move both points to locate the first one on the first ring
  p1[2] -= Dz;

  //Check bounds
  if(p1[0] >  (ringRad+detDepth)*1.001 || p1[0] < ringRad*0.999     ||
     p1[1] >  detSizeX05*1.001         || p1[1] < -detSizeX05*1.001 ||
     p1[2] >  detSizeY*1.001           || p1[2] < -1.0e-3           ){
    printf("Error: Unconsistent coincidence data.\n");
    printf(" Original P1: %15.5E %15.5E %15.5E\n",
	   p1Orig[0],
	   p1Orig[1],
	   p1Orig[2]);
    printf("Resulting P1: %15.5E %15.5E %15.5E\n",
	   p1[0], p1[1], p1[2]);
    printf("Valid transformed P1 intervals: ( [%15.5E,%15.5E], [%15.5E,%15.5E], [%15.5E,%15.5E] )\n",
	   ringRad, ringRad+detDepth,
	   -detSizeX05, detSizeX05,
	   0.0, detSizeY);
  }

  //Clamp position
  p1[0] = std::clamp(p1[0],     ringRad*1.000001, (ringRad+detDepth)*0.999999);
  p1[1] = std::clamp(p1[1], -detSizeX05*0.999999,         detSizeX05*0.999999);
  p1[2] = std::clamp(p1[2],               1.0e-6,           detSizeY*0.999999);

  //Transform the point to local detector coordinates
  p1[0] -= ringRad;
  p1[1] += detSizeX05;
  p1[2] = detSizeY - p1[2];
  
  //Return final point in local coordinates
  //Note: In local coordinates (x,y,DOI) = (p[1],p[2],p[0])
  return {p1[1], p1[2], p1[0]};
}

double attenuationFactorWithAir(const std::array<double, 3>& p1,
				const std::array<double, 3>& p2,
				double mu_cylinder,
				double mu_air,
				const std::array<double, 3>& cylinder_origin,
				double radius,
				double height) {
  
    std::array<double, 3> d = subtract(p2, p1); // Direction vector
    std::array<double, 3> o = subtract(p1, cylinder_origin); // LOR origin relative to cylinder

    // Total LOR length
    double L_total = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

    // Quadratic coefficients for intersection with infinite cylinder (aligned along Z)
    double a = d[0]*d[0] + d[1]*d[1];
    double b = 2 * (o[0]*d[0] + o[1]*d[1]);
    double c = o[0]*o[0] + o[1]*o[1] - radius*radius;

    double discriminant = b*b - 4*a*c;
    double L_cylinder = 0.0;

    if (discriminant >= 0 && a != 0.0) {
        double sqrt_disc = std::sqrt(discriminant);
        double t1 = (-b - sqrt_disc) / (2*a);
        double t2 = (-b + sqrt_disc) / (2*a);

        // Clamp to [0,1] range of LOR
        t1 = std::max(0.0, std::min(1.0, t1));
        t2 = std::max(0.0, std::min(1.0, t2));

        if (t1 != t2) {
            // Compute Z coordinates of entry and exit points
            double z1 = o[2] + d[2]*t1;
            double z2 = o[2] + d[2]*t2;

            double z_min = 0.0;
            double z_max = height;

            if (!(z1 < z_min && z2 < z_min) && !(z1 > z_max && z2 > z_max)) {
                // Clamp Z to cylinder bounds
                double tz1 = t1, tz2 = t2;
                if (z1 < z_min) tz1 = (z_min - o[2]) / d[2];
                if (z1 > z_max) tz1 = (z_max - o[2]) / d[2];
                if (z2 < z_min) tz2 = (z_min - o[2]) / d[2];
                if (z2 > z_max) tz2 = (z_max - o[2]) / d[2];

                // Final entry and exit points
                std::array<double, 3> entry = add(p1, scale(d, tz1));
                std::array<double, 3> exit  = add(p1, scale(d, tz2));

                double dx = exit[0] - entry[0];
                double dy = exit[1] - entry[1];
                double dz = exit[2] - entry[2];
                L_cylinder = std::sqrt(dx*dx + dy*dy + dz*dz);
            }
        }
    }

    double L_air = L_total - L_cylinder;
    return std::exp(mu_cylinder * L_cylinder + mu_air * L_air);
}
