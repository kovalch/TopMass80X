#include "IdeogramCombLikelihood.h"

#include "TF1.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool onlySIG  = false;
  bool useCalib = true;
  bool fSigOfMtop = false;
  
  double fSig = 0.539; // theory: 164pb-1
  //* worst case, improvable by f(w_i)
  double fCP = 2.79599875211715698e-01-0.5*p[5]; // 2.78415888547897339e-01; // prob > 0.09: 2.60629922151565552e-01;
  double fWP = 4.28881868720054626e-03+1.13365799188613892e-03+2.13942706584930420e-01-0.5*p[5]; // 6.03945180773735046e-03+1.20438553858548403e-03+2.17674896121025085e-01;
  // prob > 0.09: 6.07527932152152061e-03+1.22986722271889448e-03+2.14076578617095947e-01; // inner-branch, inter-branch, missing-jet
  double fUN = 5.01034975051879883e-01+p[5]; // 4.96665388345718384e-01; // prob > 0.09: 5.17988383769989014e-01;
  //*/

  if (onlyCP) {
    fWP = 0;
    fUN = 0;
  }
  
  if(onlySIG){
    fSig = 1.;
  }

  if(fSigOfMtop){
    double a0    = -0.00464347142282148;   // CSVT, P(chi2) weighted
    double a1    =  0.0000348026883089243; // CSVT, P(chi2) weighted
    double b0    =  477200000000; // <- Ahrens //  497700000000; // <- Kidonakis //  501000000000; // <- Langenfeld // 
    double b1    = -3322000000;  // <- Ahrens //  -3375000000;  // <- Kidonakis //  -3385000000;  // <- Langenfeld // 
    double b2    =  10050000;   // <- Ahrens //    9625000;    // <- Kidonakis //    9553000;    // <- Langenfeld // 
    double b3    = -12340;     // <- Ahrens //    -10770;     // <- Kidonakis //    -10470;     // <- Langenfeld // 
    double Lint  = 3544.844;
    double Ndata = 2418.*0.437638;

    fSig = (a0+a1*x[0])*(b0+b1*x[0]+b2*pow(x[0],2)+b3*pow(x[0],3))/pow(x[0],4)*Lint/Ndata;
  }

  //double MassOffset       =  2.70748e-01;
  //double MassSlopeMass    = -3.85914e-02;
  //double MassSlopeJES     =  2.72970e+01;
  //double MassSlopeMassJES =  0.;
  //
  //double JESOffset        =  4.84035e-03;
  //double JESSlopeMass     =  2.13560e-04;
  //double JESSlopeJES      = -2.72912e-01;
  //double JESSlopeMassJES  =  0.;

  //double MassOffset       =  5.56532e-01; // 3.51932e-01; // 
  //double MassSlopeMass    = -3.55417e-02; //-2.81632e-02; // 
  //double MassSlopeJES     =  3.17736e+01-1.27183e+01; // 2.79019e+01; // 
  //double MassSlopeMassJES = -1.12788e-01; //-1.20406e+00; // 
  //	 			    
  //double JESOffset        =  6.81159e-04; // 4.18722e-03; // 
  //double JESSlopeMass     =  3.71348e-04; // 1.93759e-04; // 
  //double JESSlopeJES      = -3.58257e-01+1.41266e-01; //-2.84423e-01; // 
  //double JESSlopeMassJES  =  1.16500e-03; // 1.41834e-02; // 
  
  if (useCalib) {

    double MassOffset0       =  1.05709e-02; //  4.61290e-01; //  4.85521e-01; // 
    double MassSlopeMass0    =  6.41444e-02; // -1.24454e-02; // -1.24608e-02; // 
    double MassSlopeJES0     = -4.61594e+00; //  2.36012e+01; //  2.35927e+01; // 
    double MassSlopeMassJES0 =  0.;          //  0.;          //  0.;          // 
			                               
    double JESOffset0        = -3.76801e-03; // -2.11106e-03; // -2.16026e-03; // 
    double JESSlopeMass0     = -1.50418e-04; //  2.39361e-04; //  2.45736e-04; // 
    double JESSlopeJES0      =  4.60147e-02; // -2.76961e-01; // -2.77323e-01; // 
    double JESSlopeMassJES0  =  0.;          //  0.;          //  0.;          // 

    double m0 = x[0];
    double j0 = x[1];

    double j1 = j0 +  JESOffset0 +  JESSlopeMass0*(m0-172.5) +  JESSlopeJES0*(j0-1.) +  JESSlopeMassJES0*(m0-172.5)*(j0-1.);
    double m1 = m0 + MassOffset0 + MassSlopeMass0*(m0-172.5) + MassSlopeJES0*(j0-1.) + MassSlopeMassJES0*(m0-172.5)*(j0-1.);

    // do additional FastSim calibration if p[6] == 1 (command line parameter -F 1)
    double MassOffset1       = (p[6] == 1) ?  1.89113 : 0.;
    double MassSlopeMass1    = 0.;
    double MassSlopeJES1     = (p[6] == 1) ? -14.2863 : 0.;
    double MassSlopeMassJES1 = 0.;
			           
    double JESOffset1        = (p[6] == 1) ? -1.98822e-02 : 0.;
    double JESSlopeMass1     = 0.;
    double JESSlopeJES1      = (p[6] == 1) ?  5.70331e-02 : 0.;
    double JESSlopeMassJES1  = 0.;
 
    x[1] = j1 +  JESOffset1 +  JESSlopeMass1*(m1-172.5) +  JESSlopeJES1*(j1-1.) +  JESSlopeMassJES1*(m1-172.5)*(j1-1.);
    x[0] = m1 + MassOffset1 + MassSlopeMass1*(m1-172.5) + MassSlopeJES1*(j1-1.) + MassSlopeMassJES1*(m1-172.5)*(j1-1.);

    //double j2 = j1 +  JESOffset1 +  JESSlopeMass1*(m1-172.5) +  JESSlopeJES1*(j1-1.) +  JESSlopeMassJES1*(m1-172.5)*(j1-1.);
    //double m2 = m1 + MassOffset1 + MassSlopeMass1*(m1-172.5) + MassSlopeJES1*(j1-1.) + MassSlopeMassJES1*(m1-172.5)*(j1-1.);
    //
    //double MassOffset2       = -2.60393e-02;
    //double MassSlopeMass2    =  2.67086e-03;
    //double MassSlopeJES2     = -2.20922e+00;
    //double MassSlopeMassJES2 =  0;
    //
    //double JESOffset2       =  1.74822e-04;
    //double JESSlopeMass2    = -3.29419e-05;
    //double JESSlopeJES2     =  2.33174e-02;
    //double JESSlopeMassJES2 =  0;
    //
    //x[1] = j2 +  JESOffset2 +  JESSlopeMass2*(m2-172.5) +  JESSlopeJES2*(j2-1.) +  JESSlopeMassJES2*(m2-172.5)*(j2-1.);
    //x[0] = m2 + MassOffset2 + MassSlopeMass2*(m2-172.5) + MassSlopeJES2*(j2-1.) + MassSlopeMassJES2*(m2-172.5)*(j2-1.);
  }
  
  //std::cout << "TEST: " << PCP (x,p) << " " << PWP (x,p) << " " << PUN (x,p) << std::endl;
  return p[0] * (fSig * (fCP * PCP (x,p) * (p[3]*PCPJES1 (x,p) + (1.-p[3])*PCPJES2 (x,p))
		       + fWP * PWP (x,p) * (p[3]*PWPJES1 (x,p) + (1.-p[3])*PWPJES2 (x,p))
		       + fUN * PUN (x,p) * (p[3]*PUNJES1 (x,p) + (1.-p[3])*PUNJES2 (x,p)))  
	    + (1.-fSig)      * PBKG(x,p) * (p[3]*PBKGJES1(x,p) + (1.-p[3])*PBKGJES2(x,p)));
  //return p[0] * (fCP * PCP(x, p) * (p[3]*PCPJES1(x, p) + (1-p[3])*PCPJES2(x,p)) + fWP * PWP(x, p) * (p[3]*PWPJES1(x, p) + (1-p[3])*PWPJES2(x,p)) + fUN * PUN(x, p) * (p[3]*PUNJES1(x, p) + (1-p[3])*PUNJES2(x, p)));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, SSVHPT
  //double q[8] = {172.62, 0.976625, 88.077, 0.796679, 7.21724, 0.0587786, 6.52763, 0.0383936};
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
  //double q[8] = {172.351, 0.994749, 85.8624, 0.819629, 7.31249, 0.0549339, 6.80791, 0.0648612};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Fall11 CSVT
  //double q[8] = {172.341, 0.996077, 85.9855, 0.814211, 7.16981, 0.0607052, 6.88698, 0.247797};
  //double q[8] = {172.321, 0.986021, 85.9632, 0.808684, 7.18716, 0.0672145, 6.85235, 0.251698};
  double q[8] = {172.324, 0.986106, 85.9343, 0.810183, 7.18309, 0.0675277, 6.85536, 0.256293}; // used for measurement
  //double q[8] = {172.328, 0.982351, 85.7756, 0.871109, 7.17985, 0.0687401, 6.68148, 0.281139}; // 5JES

  /// no probability weight
  //double q[8] = {172.063, 0.986852, 84.0195, 0.821822, 7.70989, 0.0589825, 6.92181, -0.0486285}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT


  //* hadTopMass
  double mu       = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
  //*/
  
  /* sin(theta*)
  double mu       = 1.13011e+00 - 1.15254e-02 * (x[0]-172.5) - 5.50104e-01 * (x[1]-1.);
  double sigma    = 0.3;
  
  return TMath::Landau(p[1], mu, sigma, true);
  //*/
}

//namespace cb {
//  double A(double alpha, double power) {
//    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
//  };
//  double B(double alpha, double power) {
//    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
//  };
//}
//
//double IdeogramCombLikelihood::PWP(double* x, double* p)
//{
//  //double q[12] = {172.535, 1.00874, 88.6913, -0.551851, 19.3838, 0.365133, 17.3786, 0.141755, 0.293304, 0.00829537, 0.378054, 0.0181731}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, SSVHPT
//  double q[12] = {172.167, 0.99191, 83.2081, 1.67413, 19.7897, 0.297063, 15.5488, 1.11374, 0.309901, 0.0073293, 0.441622, 0.037262}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
//
//  /// no probability weight
//  //double q[12] = {171.829, 0.876338, 86.8266, 1.37924, 19.7704, 0.262535, 17.7181, 1.15952, 0.303266, 0.00576271, 0.49986, 0.0345235}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
//
//  double N      =  1./0.01;
//  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
//  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
//  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
//  double power  =  15;
//  double t = (p[1] - mu) / sigma;
//  
//  //*
//  double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
//  double N2 = cb::A(alpha,power) / (-1.+power) * (
//                (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
//               -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
//              );
//  N = N1 + N2;
//  //*/
//  
//  if(t < alpha)
//    return 1./N * TMath::Exp(-t*t/2);
//  else
//    return 1./N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
//}
//
//double IdeogramCombLikelihood::PUN(double* x, double* p)
//{
//  //double q[12] = {175.068, 0.873118, 85.4606, 1.39263, 19.2763, 0.2894, 23.8891, 1.15432, 0.513033, 0.00981049, 0.809559, 0.0451389}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, SSVHPT
//  double q[12] = {174.484, 0.851875, 80.6582, 1.16383, 19.3731, 0.259819, 22.4076, 0.689002, 0.557451, 0.0097931, 0.805814, 0.0107206}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
//
//  /// no probability weight
//  //double q[12] = {174.85, 0.856814, 80.1727, 1.29134, 20.8871, 0.228672, 66.5038, -0.0987492, 0.667331, 0.00946654, 5.49449, 0.00248553}; // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
//
//  double N      =  1./0.01;
//  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
//  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
//  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
//  double power  =  5;
//  double t = (p[1] - mu) / sigma;
//  
//  //*
//  double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
//  double N2 = cb::A(alpha,power) / (-1.+power) * (
//                (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
//               -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
//              );
//  N = N1 + N2;
//  //*/
//  
//  if(t < alpha)
//    return 1./N * TMath::Exp(-t*t/2);
//  else
//    return 1./N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
//}

namespace MyFunctions {
  double langau(Double_t *x, Double_t *p) {
    return p[0]*(p[4]*TMath::Landau(x[0],p[1],p[2],true)+(1-p[4])*TMath::Gaus(x[0],p[1],p[3],true));
  }
}

double IdeogramCombLikelihood::PWP(double* x, double* p)
{
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
  //double q[12] = {177.506, 0.72069, 72.7239, -0.417889, 23.0854, 0.141175, 4.41442, 0.951157, 8.21835, -0.00774364, -2.46418, -0.665209}; 
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
  //double q[12] = {177.506, 0.72069, 72.7239, 0., 23.0854, 0.141175, 4.41442, 0., 8.21835, -0.00774364, -2.46418, 0.};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Summer11 178, 181, CSVT
  //double q[12] = {175.998, 0.680979, 72.8323, 0.0836664, 21.5934, 0.0731161, 6.49105, 0.642425, 7.49304, -0.0451434, -1.27007, 0.300344};
  //double q[12] = {175.917, 0.662955, 73.7242, 0.361619, 21.6296, 0.103931, 6.6901, 0.537813, 7.45372, -0.0525605, -1.74566, 0.120026};
  double q[12] = {175.941, 0.661709, 73.1694, 0.36766, 21.61, 0.100468, 6.74084, 0.561399, 7.44604, -0.0502921, -1.84577, 0.100049}; // used for measurement
  //double q[12] = {175.979, 0.646089, 73.2551, 0.322292, 21.6434, 0.104014, 5.95933, 0.611622, 7.52061, -0.0627159, -2.45757, 0.496283}; // 5JES

  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = 0.9;

  TF1 func = TF1("func", MyFunctions::langau, 100., 550., 5);
  
  func.SetParameters(1., mu, sigmaL, sigmaG, ratio);
    
  ScanPoint currentPoint = std::make_pair(x[0],x[1]);
  ScanPointMap::const_iterator iter = PWPintegrals_.find(currentPoint);
  double currentIntegral = -1.;
  if(iter == PWPintegrals_.end()){
    currentIntegral = func.Integral(100.,550.);
    PWPintegrals_[currentPoint] = currentIntegral;
  }
  else{
    currentIntegral = iter->second;
  }
  double N = 1. / currentIntegral; // func.Integral(100.,550.); // 

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihood::PUN(double* x, double* p)
{
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
  //double q[12] = {176.206, 0.697207, 69.1265, 1.88218, 22.2574, 0.18319, 16.9791, 0.644276, 9.61471, -0.0348847, 1.59068, 1.18297};
  // prob > 0.09 & dRbb > 1.5, FALL11 161, 172, 184, CSVT
  //double q[12] = {176.206, 0.697207, 69.1265, 0., 22.2574, 0.18319, 16.9791, 0., 9.61471, -0.0348847, 1.59068, 0.};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Summer11 178, 181, CSVT
  //double q[12] = {175.43, 0.756711, 70.1407, 1.33444, 21.1268, 0.197732, 15.356, 1.36539, 8.85433, 0.00154935, 0.188517, 0.757486};
  //double q[12] = {175.51, 0.731126, 70.2508, 1.80935, 21.2712, 0.182935, 15.1928, 1.71509, 8.95565, -0.0277434, -0.739389, 0.751678};
  double q[12] = {175.453, 0.736703, 70.2111, 1.69023, 21.1899, 0.189805, 15.3274, 1.60513, 8.87531, -0.0238031, -0.117816, 0.760937}; // used for measurement
  //double q[12] = {175.441, 0.736269, 69.6938, 1.70864, 21.1644, 0.194426, 14.8559, 1.75101, 8.82208, -0.0224898, -0.222645, 0.658529}; // 5JES

  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = 0.8;

  TF1 func = TF1("func", MyFunctions::langau, 100., 550., 5);
  
  func.SetParameters(1., mu, sigmaL, sigmaG, ratio);
    
  ScanPoint currentPoint = std::make_pair(x[0],x[1]);
  ScanPointMap::const_iterator iter = PUNintegrals_.find(currentPoint);
  double currentIntegral = -1.;
  if(iter == PUNintegrals_.end()){
    currentIntegral = func.Integral(100.,550.);
    PUNintegrals_[currentPoint] = currentIntegral;
  }
  else{
    currentIntegral = iter->second;
  }
  double N = 1. / currentIntegral; // func.Integral(100.,550.); // 

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihood::PCPJES1(double *x, double *p)
{
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, SSVHPT
  //double q[12] = {84.5383, -0.01724, 58.0165, -0.0536437, 4.39805, -0.00844352, 11.6747, -0.0682757, 5.5907, 0.00266851, -6.33346, -0.0166725};
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT
  //double q[12] = {84.4836, -0.00235611, 57.7472, 0.260235, 4.3087, -0.00226071, 12.085, 0.108488, 5.43718, -0.0112025, -7.38103, -0.125601};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Summer11 178, 181, CSVT
  //double q[12] = {84.4708, -0.00114955, 57.3044, 0.354774, 4.25771, 0.00114764, 11.86, 0.18864, 5.42256, -0.00984776, -7.03663, -0.219064};
  //double q[12] = {84.4567, -0.00512651, 57.1306, 0.179133, 4.24987, 0.000300746, 11.7124, 0.0582619, 5.44599, -0.00559353, -6.87845, -0.0609114};
  double q[12] = {84.4587, -0.00524409, 57.1528, 0.178028, 4.25284, 0.000167644, 11.7426, 0.0582509, 5.43288, -0.00549247, -6.9163, -0.0576681}; // used for measurement
  //double q[12] = {84.4714, -0.00288904, 57.2879, 0.221738, 4.25403, 0.00163248, 11.8082, 0.0886389, 5.43406, -0.00640861, -6.93172, -0.0742648}; // 5JES

  /// no probability weight
  //double q[12] = {84.5916, -0.00171283, 91.7268, -0.0664178, 5.08861, -0.00321848, 24.5077, -0.0313107, 6.48301, -0.0100923, -20.077, 0.00596117}; // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT

  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PWPJES1(double *x, double *p)
{
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, SSVHPT
  //double q[12] = {84.0048, 0.0204403, 27.3415, -0.959433, 4.34126, 0.00500131, 5.28954, -0.776945, 6.50917, -0.0179892, -6.64902, 0.53598};
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT
  //double q[12] = {84.0586, 0.0314863, 28.6706, -0.524549, 4.31361, 0.00693565, 5.86459, -0.361966, 6.27985, -0.0230866, -6.57438, 0.40124};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Summer11 178, 181, CSVT
  //double q[12] = {84.0542, 0.0199065, 27.7822, -0.0726202, 4.2925, -0.0032254, 5.41803, -0.130014, 6.27911, -0.0178925, -6.36262, 0.0745195};
  //double q[12] = {84.0441, 0.0102286, 27.5264, -0.133357, 4.2872, -0.00758038, 5.24515, -0.216303, 6.29317, -0.0168775, -6.20706, 0.163248};
  double q[12] = {84.042, 0.0111381, 27.5838, -0.179698, 4.28814, -0.00702476, 5.27829, -0.234302, 6.28184, -0.0166415, -6.23916, 0.181574}; // used for measurement
  //double q[12] = {84.059, 0.00363846, 27.6651, -0.0122909, 4.28601, -0.0106779, 5.2756, -0.18915, 6.2658, -0.0133067, -6.2384, 0.0742743}; // 5JES

  /// no probability weight
  //double q[12] = {84.8324, 0.0162648, 48.9022, 0.274153, 5.32539, -0.00247155, 13.0281, 0.0295863, 7.5104, -0.0199612, -17.7564, -0.177945}; // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT
  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PUNJES1(double *x, double *p)
{
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, SSVHPT
  //double q[12] = {84.7265, 0.00352151, 19.5715, -0.0682227, 4.85665, 0.00148563, 3.79153, -0.0669777, 7.24342, 0.00743381, -8.78048, 0.0662842};
  // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT
  //double q[12] = {84.7383, 0.00308977, 16.0727, 0.38017, 4.84492, 0.00678078, 2.43329, 0.14005, 7.04282, 0.00211119, -6.0521, -0.17391};
  // prob > 0.09 & dRbb > 1.5, ROOFIT, Summer11 178, 181, CSVT
  //double q[12] = {84.634, 0.0265456, 16.4579, 0.038807, 4.77021, 0.0176135, 2.39677, 0.0458818, 7.08442, -0.0108841, -6.42561, 0.0666976};
  //double q[12] = {84.6112, 0.017789, 15.7834, -0.181863, 4.76847, 0.0184747, 2.04143, -0.105465, 7.14669, -0.00782099, -6.08471, 0.132059};
  double q[12] = {84.6224, 0.0203645, 16.0195, -0.293667, 4.77153, 0.0196543, 2.16641, -0.150288, 7.09421, -0.00659033, -6.2366, 0.235723}; // used for measurement
  //double q[12] = {84.6333, 0.0205799, 16.0275, -0.280052, 4.76909, 0.018934, 2.18012, -0.12839, 7.08441, -0.00812571, -6.2732, 0.22637}; // 5JES

  /// no probability weight
  //double q[12] = {86.7914, 0.00700164, 26.5565, 0.698038, 6.48728, 0.00675855, 5.76356, 0.328417, 7.74975, -0.00159438, -13.9122, -0.274122}; // prob > 0.09 & dRbb > 1.5 - mean top mass, FALL11 161, 172, 184, CSVT
  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PCPJES2(double *x, double *p)
{
  double q[12] = {83.2889, 0.0136309, 60.2892, 0.750832, 5.60129, 0.00022743, 11.0848, 0.522803, 8.13666, -0.0095548, -6.56246, -0.184139}; // prob > 0.01 & dRbb > 1.5
  //double q[12] = {83.2947, 0.012297, 59.775, 0.523396, 5.49951, 0.00162267, 11.7272, 0.391729, 7.96111, -0.0104956, -7.21539, -0.165301}; // prob > 0.09
  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PWPJES2(double *x, double *p)
{
  double q[12] = {82.6176, -0.0343541, 39.8918, -1.25747, 5.69796, -0.0139934, 11.665, -0.815988, 9.27449, 0.0224656, -11.8063, 0.211439}; // prob > 0.01 & dRbb > 1.5
  //double q[12] = {82.6203, -0.0444749, 37.7678, -0.806399, 5.56005, -0.0171444, 10.5941, -0.647278, 9.05804, 0.0314525, -11.9595, 0.0656524}; // prob > 0.09
  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PUNJES2(double *x, double *p)
{
  double q[12] = {82.7644, -0.0131422, 25.0316, 0.809681, 6.32278, -0.0197688, 5.99963, 0.479629, 11.361, 0.00848729, -12.3799, -0.414559}; // prob > 0.01 & dRbb > 1.5
  //double q[12] = {82.9722, -0.0108677, 22.7339, 1.22133, 6.21211, -0.021831, 4.63532, 0.681988, 10.6728, 0.00302816, -11.5825, -0.562818}; // prob > 0.09
  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PJES(double* x, double* p, double *q)
{
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
   
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  static const double sqrt2 = sqrt(2);
  double N1     =  sqrt2*sigma1;
  double N2     =  sqrt2*sigma2;
  
  double N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
}

namespace MyFunctions {
  double gammaland(Double_t *x, Double_t *p) {
    if(x[0] < p[2]) return 0.;
    else return (p[0]*TMath::GammaDist(x[0],p[1],p[2],p[3])+p[4]*TMath::Landau(x[0],p[5],p[6]))/9.844;
  }
}
double IdeogramCombLikelihood::PBKG(double *x, double *p)
{
  // [0]*TMath::GammaDist(x,[1],[2],[3])+[4]*TMath::Landau(x,[5],[6])
  //double q[7] = {8.07678, 4.06005, 89.9024, 33.6036, 0.0881897, 211.786, 24.4228}; // prob > 0.09 & dRbb > 1.5, ???TCHEM???
  //double q[7] = {8.56530, 4.00234, 90.5065, 33.3625, 0.0564579, 221.309, 29.3280}; // prob > 0.09 & dRbb > 1.5, CSVT
  //double q[7] = {6.54902, 3.63210, 126.189 , 32.2675, 0.168945 , 158.400, 22.1696}; // prob > 0.09 & dRbb > 1.5, CSVT
  double q[7] = {8.37699, 4.30272,  87.9676, 32.6083, 0.0628411, 220.731, 29.9386}; // prob > 0.09 & dRbb > 1.5, CSVT
  //double q[7] = {8.67304, 4.80112, 87.5538, 28.0572, 0.0689699, 215.102, 21.1638}; // OBTAINED FROM REMIXING SIGNAL MC ONLY
  //double q[7] = {8.19217, 5.16177, 80.6748, 28.0975, 0.0795278, 220.641, 25.6501}; // OBTAINED WITH MET CUT

  TF1 func = TF1("func", MyFunctions::gammaland, 100., 550., 7);
  
  func.SetParameters(q);
  func.SetParameter(1, p[4]*q[1]);
    
  if(PBKGintegral_ < 0.){
    PBKGintegral_ = func.Integral(100.,550.);
  }
  double N = 1. / PBKGintegral_; // func.Integral(100.,550.); // 

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihood::PBKGJES1(double *x, double *p)
{
  //double q[3] = {88.1473, 6.3534, 7.8005}; // prob > 0.09 & dRbb > 1.5 - mean top mass, ???TCHEM???
  //double q[3] = {85.2035, 7.08030, 12.0783}; // prob > 0.09, ???TCHEM???
  //double q[3] = {86.8970, 5.79471, 7.22339}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  //double q[3] = {87.1875, 6.29255, 7.36406}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  double q[3] = {87.1042, 5.83934, 7.02516}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  //double q[3] = {86.2543, 5.59517, 7.05928}; // OBTAINED FROM REMIXING SIGNAL MC ONLY
  //double q[3] = {87.0748, 5.95375, 6.75081}; // OBTAINED WITH MET CUT

  return PBKGJES(x, p, q);
}

double IdeogramCombLikelihood::PBKGJES2(double *x, double *p)
{
  double q[3] = {85.4051, 7.02919, 11.9709};
  return PBKGJES(x, p, q);
}

double IdeogramCombLikelihood::PBKGJES(double *x, double *p, double *q)
{
  double mu     = q[0];
  double sigma1 = q[1];
  double sigma2 = q[2];
   
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  static const double sqrt2 = sqrt(2);
  double N1     =  sqrt2*sigma1;
  double N2     =  sqrt2*sigma2;
  
  double N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
}
