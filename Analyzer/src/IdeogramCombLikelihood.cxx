#include "IdeogramCombLikelihood.h"

//#include <iostream>

#include "TF1.h"
#include "TMath.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool onlySIG  = false;
  bool useCalib = true;
  bool fSigOfMtop = false;
  
  double fSig = 0.539; // theory: 164pb-1
  //* worst case, improvable by f(w_i)
  double fCP = 2.31859251856803894e-01-0.5*p[5]; //2.79599875211715698e-01-0.5*p[5];
  double fWP = 4.23728721216320992e-03+9.87463630735874176e-04+2.18985587358474731e-01-0.5*p[5]; //4.28881868720054626e-03+1.13365799188613892e-03+2.13942706584930420e-01-0.5*p[5];
  // inner-branch, inter-branch, missing-jet
  double fUN = 5.43930411338806152e-01+p[5]; //5.01034975051879883e-01+p[5];
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

  if (useCalib) {

    double MassOffset0       = 0.;
    double MassSlopeMass0    = 0.;
    double MassSlopeJES0     = 0.;
    double MassSlopeMassJES0 = 0.;
			                               
    double JESOffset0        = 0.;
    double JESSlopeMass0     = 0.;
    double JESSlopeJES0      = 0.;
    double JESSlopeMassJES0  = 0.;

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
  
  return p[0] * (fSig * (fCP * PCP (x,p) * (p[3]*PCPJES1 (x,p) + (1.-p[3])*PCPJES2 (x,p))
		               + fWP * PWP (x,p) * (p[3]*PWPJES1 (x,p) + (1.-p[3])*PWPJES2 (x,p))
		               + fUN * PUN (x,p) * (p[3]*PUNJES1 (x,p) + (1.-p[3])*PUNJES2 (x,p)))
	       + (1.-fSig)       * PBKG(x,p) * (p[3]*PBKGJES1(x,p) + (1.-p[3])*PBKGJES2(x,p)));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  //double q[8] = {172.324, 0.986106, 85.9343, 0.810183, 7.18309, 0.0675277, 6.85536, 0.256293}; // used for measurement
  double q[8] = {171.316, 0.964353, 82.0417, 0.669057, 7.17516, 0.0635755, 5.71509, 0.0795821}; // L7Correction

  double mu       = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
 }

namespace MyFunctions {
  double langau(Double_t *x, Double_t *p) {
    return p[0]*(p[4]*TMath::Landau(x[0],p[1],p[2],true)+(1-p[4])*TMath::Gaus(x[0],p[1],p[3],true));
  }
  double lanvog(Double_t *x, Double_t *p) {
    return p[0]*(p[5]*TMath::Landau(x[0],p[1],p[3],true)+(1-p[5])*TMath::Gaus(x[0],p[2],p[4],true));
  }
}

double IdeogramCombLikelihood::PWP(double* x, double* p)
{
  //double q[12] = {175.941, 0.661709, 73.1694, 0.36766, 21.61, 0.100468, 6.74084, 0.561399, 7.44604, -0.0502921, -1.84577, 0.100049}; // used for measurement
  double q[16] = {176.044, 0.249872, 84.6797, 0.0325573, 171.527, 0.98991, 87.0197, 0.0133321, 22.2626, -0.0108761, 16.1794, 0.398061, 8.92348, 0.0677935, 11.293, -2.7555}; // L7Correction

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = 8.11788e-01;

  //TF1 func = TF1("func", MyFunctions::langau, 100., 550., 5);
  TF1 func = TF1("func", MyFunctions::lanvog, 100., 550., 6);
  
  func.SetParameters(1., muL, muG, sigmaL, sigmaG, ratio);
    
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
  //double q[12] = {175.453, 0.736703, 70.2111, 1.69023, 21.1899, 0.189805, 15.3274, 1.60513, 8.87531, -0.0238031, -0.117816, 0.760937}; // used for measurement
  double q[16] = {175.133, 0.271484, 35.1462, 0.774396, 173.809, 0.946292, 77.146, 1.65189, 21.0671, 0.0181014, 0.428224, 0.666531, 9.45852, 0.0400708, -11.9838, 0.0918909}; // L7Correction

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = 7.72757e-01;

  //TF1 func = TF1("func", MyFunctions::langau, 100., 550., 5);
  TF1 func = TF1("func", MyFunctions::lanvog, 100., 550., 6);

  func.SetParameters(1., muL, muG, sigmaL, sigmaG, ratio);
    
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
  //double q[12] = {84.4587, -0.00524409, 57.1528, 0.178028, 4.25284, 0.000167644, 11.7426, 0.0582509, 5.43288, -0.00549247, -6.9163, -0.0576681}; // used for measurement
  //double q[12] = {80.7147, -0.0188437, 57.4631, 0.260386, 3.5588, -0.00823717, 11.1024, 0.0370445, 5.81416, 0.00514867, -8.59795, -0.279046}; // L7Correction
  double q[12] = {86.0752, -0.0192945, 61.9831, 0.23695, 3.92178, -0.00866597, 11.9513, 0.0231462, 6.4164, 0.00638338, -8.47745, -0.252667};

  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PWPJES1(double *x, double *p)
{
  //double q[12] = {84.042, 0.0111381, 27.5838, -0.179698, 4.28814, -0.00702476, 5.27829, -0.234302, 6.28184, -0.0166415, -6.23916, 0.181574}; // used for measurement
  //double q[12] = {82.2216, -0.00209122, 22.518, -0.191243, 4.08095, -0.00377955, 0.464987, -0.205658, 6.7557, -0.00202161, -10.7312, 0.0934139}; // L7Correction
  double q[12] = {87.4292, 0.0040691, 24.0689, 0.0386941, 4.41083, -0.00307802, 0.962317, -0.0797275, 7.3466, -0.00833658, -10.0061, -0.0529531};

  return PJES(x, p, q);
}

double IdeogramCombLikelihood::PUNJES1(double *x, double *p)
{
  //double q[12] = {84.6224, 0.0203645, 16.0195, -0.293667, 4.77153, 0.0196543, 2.16641, -0.150288, 7.09421, -0.00659033, -6.2366, 0.235723}; // used for measurement
  //double q[12] = {83.407, 0.0197121, 11.6287, 0.157607, 4.6822, 0.0100276, -0.205017, 0.23718, 7.50906, -0.0111461, -7.79563, 0.00502378}; // L7Correction
  double q[12] = {88.502, 0.0257782, 12.5532, 0.0892862, 4.97979, 0.0125393, 0.108106, 0.220309, 8.28809, -0.0155682, -7.42246, 0.0794568};

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
  //double q[7] = {8.37699, 4.30272,  87.9676, 32.6083, 0.0628411, 220.731, 29.9386}; // prob > 0.09 & dRbb > 1.5, CSVT
  double q[7] = {8.37699/*0.852647*/, 3.88501, 93.2179, 33.7309, 0.0628411/*0.147353*/, 206.862, 25.7198}; // L7Correction

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
  //double q[3] = {87.1042, 5.83934, 7.02516}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  double q[3] = {91.9298, 6.23531, 8.12332}; // L7Correction

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
