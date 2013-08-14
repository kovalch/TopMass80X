#include "IdeogramCombLikelihoodAllJets.h"

#include "Helper.h"
#include "ProgramOptionsReader.h"

//#include <iostream>

#include "TF1.h"
#include "TFile.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

typedef ProgramOptionsReader po;

IdeogramCombLikelihoodAllJets::IdeogramCombLikelihoodAllJets() :
    PBKGintegral_(-1)
{
  int channelID = Helper::channelID();

  TString fileName = "Calibration_";
  if     (channelID == Helper::kAllJets) fileName += "AllJets";
  else if(channelID == Helper::kLeptonJets || channelID == Helper::kMuonJets || channelID == Helper::kElectronJets) fileName += "LeptonJets";
  fileName += ".root";
  TFile* calibrationFile = TFile::Open(fileName);
  RooWorkspace* cw = (RooWorkspace*)calibrationFile->Get("workspaceMtop")->Clone();

  parsCP_     = new double[8];
  parsWP_     = new double[17];
  parsUN_     = new double[17];
  parsCPJES_  = new double[12];
  parsWPJES_  = new double[12];
  parsUNJES_  = new double[12];
  parsBKG_    = new double[8];
  parsBKGJES_ = new double[3];

  // parameters for mTop correct permutations
  parsCP_[0] = cw->var("par_0_0")->getVal(); parsCP_[1] = cw->var("par_0_1")->getVal(); parsCP_[2] = cw->var("par_0_2")->getVal(); parsCP_[3] = cw->var("par_0_3")->getVal();
  parsCP_[4] = cw->var("par_0_4")->getVal(); parsCP_[5] = cw->var("par_0_5")->getVal(); parsCP_[6] = cw->var("par_0_6")->getVal(); parsCP_[7] = cw->var("par_0_7")->getVal();

  // parameters for mTop wrong parmutations
  parsWP_[ 0] = cw->var("par_1_0")->getVal();  parsWP_[ 1] = cw->var("par_1_1")->getVal();  parsWP_[ 2] = cw->var("par_1_2")->getVal();  parsWP_[ 3] = cw->var("par_1_3")->getVal();
  parsWP_[ 4] = cw->var("par_1_4")->getVal();  parsWP_[ 5] = cw->var("par_1_5")->getVal();  parsWP_[ 6] = cw->var("par_1_6")->getVal();  parsWP_[ 7] = cw->var("par_1_7")->getVal();
  parsWP_[ 8] = cw->var("par_1_8")->getVal();  parsWP_[ 9] = cw->var("par_1_9")->getVal();  parsWP_[10] = cw->var("par_1_10")->getVal(); parsWP_[11] = cw->var("par_1_11")->getVal();
  parsWP_[12] = cw->var("par_1_12")->getVal(); parsWP_[13] = cw->var("par_1_13")->getVal(); parsWP_[14] = cw->var("par_1_14")->getVal(); parsWP_[15] = cw->var("par_1_15")->getVal();
  parsWP_[16] = cw->var("ratio_1")->getVal();

  // parameters for mTop unmatched permutations
  parsUN_[ 0] = cw->var("par_2_0")->getVal();  parsUN_[ 1] = cw->var("par_2_1")->getVal();  parsUN_[ 2] = cw->var("par_2_2")->getVal();  parsUN_[ 3] = cw->var("par_2_3")->getVal();
  parsUN_[ 4] = cw->var("par_2_4")->getVal();  parsUN_[ 5] = cw->var("par_2_5")->getVal();  parsUN_[ 6] = cw->var("par_2_6")->getVal();  parsUN_[ 7] = cw->var("par_2_7")->getVal();
  parsUN_[ 8] = cw->var("par_2_8")->getVal();  parsUN_[ 9] = cw->var("par_2_9")->getVal();  parsUN_[10] = cw->var("par_2_10")->getVal(); parsUN_[11] = cw->var("par_2_11")->getVal();
  parsUN_[12] = cw->var("par_2_12")->getVal(); parsUN_[13] = cw->var("par_2_13")->getVal(); parsUN_[14] = cw->var("par_2_14")->getVal(); parsUN_[15] = cw->var("par_2_15")->getVal();
  parsUN_[16] = cw->var("ratio_2")->getVal();

  // parameters for JES correct permutations
  parsCPJES_[0] = cw->var("par_3_0")->getVal(); parsCPJES_[1] = cw->var("par_3_1")->getVal(); parsCPJES_[ 2] = cw->var("par_3_2")->getVal();  parsCPJES_[ 3] = cw->var("par_3_3")->getVal();
  parsCPJES_[4] = cw->var("par_3_4")->getVal(); parsCPJES_[5] = cw->var("par_3_5")->getVal(); parsCPJES_[ 6] = cw->var("par_3_6")->getVal();  parsCPJES_[ 7] = cw->var("par_3_7")->getVal();
  parsCPJES_[8] = cw->var("par_3_8")->getVal(); parsCPJES_[9] = cw->var("par_3_9")->getVal(); parsCPJES_[10] = cw->var("par_3_10")->getVal(); parsCPJES_[11] = cw->var("par_3_11")->getVal();

  // parameters for JES wrong permutations
  parsWPJES_[0] = cw->var("par_4_0")->getVal(); parsWPJES_[1] = cw->var("par_4_1")->getVal(); parsWPJES_[ 2] = cw->var("par_4_2")->getVal();  parsWPJES_[ 3] = cw->var("par_4_3")->getVal();
  parsWPJES_[4] = cw->var("par_4_4")->getVal(); parsWPJES_[5] = cw->var("par_4_5")->getVal(); parsWPJES_[ 6] = cw->var("par_4_6")->getVal();  parsWPJES_[ 7] = cw->var("par_4_7")->getVal();
  parsWPJES_[8] = cw->var("par_4_8")->getVal(); parsWPJES_[9] = cw->var("par_4_9")->getVal(); parsWPJES_[10] = cw->var("par_4_10")->getVal(); parsWPJES_[11] = cw->var("par_4_11")->getVal();

  // parameters for JES unmatched permutations
  parsUNJES_[0] = cw->var("par_5_0")->getVal(); parsUNJES_[1] = cw->var("par_5_1")->getVal(); parsUNJES_[ 2] = cw->var("par_5_2")->getVal();  parsUNJES_[ 3] = cw->var("par_5_3")->getVal();
  parsUNJES_[4] = cw->var("par_5_4")->getVal(); parsUNJES_[5] = cw->var("par_5_5")->getVal(); parsUNJES_[ 6] = cw->var("par_5_6")->getVal();  parsUNJES_[ 7] = cw->var("par_5_7")->getVal();
  parsUNJES_[8] = cw->var("par_5_8")->getVal(); parsUNJES_[9] = cw->var("par_5_9")->getVal(); parsUNJES_[10] = cw->var("par_5_10")->getVal(); parsUNJES_[11] = cw->var("par_5_11")->getVal();

  // parameters for mTop background
  parsBKG_[0] = cw->var("topBKGGammaNorm")->getVal();
  parsBKG_[1] = cw->var("topBKGGammaGamma")->getVal(); parsBKG_[2] = cw->var("topBKGGammaMu")->getVal(); parsBKG_[3] = cw->var("topBKGGammaBeta")->getVal();
  parsBKG_[4] = cw->var("topBKGLandauMean")->getVal(); parsBKG_[5] = cw->var("topBKGLandauSigma")->getVal();
  parsBKG_[6] = -1; parsBKG_[7] = -1;

  // parameters for JES background
  parsBKGJES_[0] = cw->var("wBKGMean")->getVal(); parsBKGJES_[1] = cw->var("wBKGSigmaLeft")->getVal(); parsBKGJES_[2] = cw->var("wBKGSigmaRight")->getVal();

  fSig_ = cw->var("fSig")->getVal(); // 0.630491408; // theory: 164pb-1
  fCP_  = cw->var("fCP")->getVal();  // 0.233318433; //
  fWP_  = cw->var("fWP")->getVal();  // 0.16536065;  //
  fUN_  = cw->var("fUN")->getVal();  // 0.601320933; //

  calibrationFile->Close();
}

double IdeogramCombLikelihoodAllJets::Evaluate(double *x, double *p) {
  bool useCalib   = true;
  bool onlyCP     = false;
  bool onlySIG    = false;
  bool fSigOfMtop = false;
  
  fCP_ -= 0.5*p[5];
  fWP_ -= 0.5*p[5];
  fUN_ += p[5];

  if (onlyCP) {
    fWP_ = 0;
    fUN_ = 0;
  }
  
  if(onlySIG){
    fSig_ = 1.;
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

    fSig_ = (a0+a1*x[0])*(b0+b1*x[0]+b2*pow(x[0],2)+b3*pow(x[0],3))/pow(x[0],4)*Lint/Ndata;
  }

  if (useCalib) {

    double MassOffset0       = 0.; // 3.31715e-01;
    double MassSlopeMass0    = 0.; // 8.22606e-03;
    double MassSlopeJES0     = 0.; // 4.60359e-01;
    double MassSlopeMassJES0 = 0.; //-2.22772e-01;

    double JESOffset0        = 0.; //-8.78984e-04;
    double JESSlopeMass0     = 0.; //-4.66983e-05;
    double JESSlopeJES0      = 0.; //-5.75365e-02;
    double JESSlopeMassJES0  = 0.; // 1.56276e-03;

    double m0 = x[0];
    double j0 = x[1];

    double j1 = j0 +  JESOffset0 +  JESSlopeMass0*(m0-172.5) +  JESSlopeJES0*(j0-1.) +  JESSlopeMassJES0*(m0-172.5)*(j0-1.);
    double m1 = m0 + MassOffset0 + MassSlopeMass0*(m0-172.5) + MassSlopeJES0*(j0-1.) + MassSlopeMassJES0*(m0-172.5)*(j0-1.);

    // do additional FastSim calibration if p[6] == 1 (command line parameter -F 1)
    double MassOffset1       = (p[6] == 1) ?  1.83190 : 0.;
    double MassSlopeMass1    = 0.;
    double MassSlopeJES1     = (p[6] == 1) ? -5.36113 : 0.;
    double MassSlopeMassJES1 = 0.;
			           
    double JESOffset1        = (p[6] == 1) ? -2.11957e-02 : 0.;
    double JESSlopeMass1     = 0.;
    double JESSlopeJES1      = (p[6] == 1) ?  3.15759e-02 : 0.;
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
  
  return p[0] * (fSig_ * (fCP_ * PCP (x,p) * (p[3]*PCPJES1 (x,p))   // + (1.-p[3])*PCPJES2 (x,p))
		               +  fWP_ * PWP (x,p) * (p[3]*PWPJES1 (x,p))   // + (1.-p[3])*PWPJES2 (x,p))
		               +  fUN_ * PUN (x,p) * (p[3]*PUNJES1 (x,p)))  // + (1.-p[3])*PUNJES2 (x,p)))
	       + (1.-fSig_)        * PBKG(x,p) * (p[3]*PBKGJES1(x,p))); // + (1.-p[3])*PBKGJES2(x,p)));
}


double IdeogramCombLikelihoodAllJets::PCP(double *x, double *p) {
  //double q[8] = {172.324, 0.986106, 85.9343, 0.810183, 7.18309, 0.0675277, 6.85536, 0.256293}; // used for measurement
  //double q[8] = {171.321, 0.96406, 82.4588, 0.667791, 7.17632, 0.0636133, 5.65153, 0.0753321}; // L7Correction
  double *q = parsCP_;

  double mu    = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
 }

double IdeogramCombLikelihoodAllJets::PWP(double* x, double* p)
{
  //double q[12] = {175.941, 0.661709, 73.1694, 0.36766, 21.61, 0.100468, 6.74084, 0.561399, 7.44604, -0.0502921, -1.84577, 0.100049}; // used for measurement
  //double q[16] = {166.62, 0.479989, 65.8855, 0.817456, 172.708, 1.0786, 92.6114, -1.10235, 15.9693, -0.0223771, 0.348473, 0.177825, 7.08448, 0.115902, 6.08694, -1.16462}; // L7Correction
  double *q = parsWP_;

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = q[16]; //8.11788e-01;

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
  double N = 1. / currentIntegral;

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihoodAllJets::PUN(double* x, double* p)
{
  //double q[12] = {175.453, 0.736703, 70.2111, 1.69023, 21.1899, 0.189805, 15.3274, 1.60513, 8.87531, -0.0238031, -0.117816, 0.760937}; // used for measurement
  //double q[16] = {178.286, 0.26325, 45.3131, 0.692725, 173.347, 0.951113, 79.5679, 1.57497, 22.5683, 0.0117107, 7.02635, 0.604059, 9.93498, 0.0310997, -12.2806, 0.000101749}; // L7Correction
  double *q = parsUN_;

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = q[16]; //7.72757e-01;

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
  double N = 1. / currentIntegral;

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihoodAllJets::PCPJES1(double *x, double *p)
{
  //double q[12] = {84.4587, -0.00524409, 57.1528, 0.178028, 4.25284, 0.000167644, 11.7426, 0.0582509, 5.43288, -0.00549247, -6.9163, -0.0576681}; // used for measurement
  //double q[12] = {80.7172, -0.0192033, 58.9908, 0.216647, 3.56043, -0.00835824, 11.4794, 0.0100302, 5.81594, 0.00560078, -8.77514, -0.240875}; // L7Correction
  double *q = parsCPJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PWPJES1(double *x, double *p)
{
  //double q[12] = {84.042, 0.0111381, 27.5838, -0.179698, 4.28814, -0.00702476, 5.27829, -0.234302, 6.28184, -0.0166415, -6.23916, 0.181574}; // used for measurement
  //double q[12] = {82.6125, 0.0113286, 14.1711, 0.94118, 4.20692, 0.000889332, -1.78333, 0.300516, 6.74883, -0.0193121, -6.19779, -0.343831}; // L7Correction
  double *q = parsWPJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PUNJES1(double *x, double *p)
{
  //double q[12] = {84.6224, 0.0203645, 16.0195, -0.293667, 4.77153, 0.0196543, 2.16641, -0.150288, 7.09421, -0.00659033, -6.2366, 0.235723}; // used for measurement
  //double q[12] = {83.069, 0.0126392, 15.8747, 0.101566, 4.56234, 0.00662185, 0.145635, 0.188651, 7.54043, -0.00766768, -10.0228, 0.00347318}; // L7Correction
  double *q = parsUNJES_;

  return PJES(x, p, q);
}

//double IdeogramCombLikelihoodAllJets::PCPJES2(double *x, double *p)
//{
//  double q[12] = {83.2889, 0.0136309, 60.2892, 0.750832, 5.60129, 0.00022743, 11.0848, 0.522803, 8.13666, -0.0095548, -6.56246, -0.184139}; // prob > 0.01 & dRbb > 1.5
//  //double q[12] = {83.2947, 0.012297, 59.775, 0.523396, 5.49951, 0.00162267, 11.7272, 0.391729, 7.96111, -0.0104956, -7.21539, -0.165301}; // prob > 0.09
//  return PJES(x, p, q);
//}
//
//double IdeogramCombLikelihoodAllJets::PWPJES2(double *x, double *p)
//{
//  double q[12] = {82.6176, -0.0343541, 39.8918, -1.25747, 5.69796, -0.0139934, 11.665, -0.815988, 9.27449, 0.0224656, -11.8063, 0.211439}; // prob > 0.01 & dRbb > 1.5
//  //double q[12] = {82.6203, -0.0444749, 37.7678, -0.806399, 5.56005, -0.0171444, 10.5941, -0.647278, 9.05804, 0.0314525, -11.9595, 0.0656524}; // prob > 0.09
//  return PJES(x, p, q);
//}
//
//double IdeogramCombLikelihoodAllJets::PUNJES2(double *x, double *p)
//{
//  double q[12] = {82.7644, -0.0131422, 25.0316, 0.809681, 6.32278, -0.0197688, 5.99963, 0.479629, 11.361, 0.00848729, -12.3799, -0.414559}; // prob > 0.01 & dRbb > 1.5
//  //double q[12] = {82.9722, -0.0108677, 22.7339, 1.22133, 6.21211, -0.021831, 4.63532, 0.681988, 10.6728, 0.00302816, -11.5825, -0.562818}; // prob > 0.09
//  return PJES(x, p, q);
//}

double IdeogramCombLikelihoodAllJets::PJES(double* x, double* p, double *q)
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

double IdeogramCombLikelihoodAllJets::PBKG(double *x, double *p)
{
  // [0]*TMath::GammaDist(x,[1],[2],[3])+[4]*TMath::Landau(x,[5],[6])
  //double q[7] = {8.37699, 4.30272,  87.9676, 32.6083, 0.0628411, 220.731, 29.9386}; // prob > 0.09 & dRbb > 1.5, CSVT
  //double q[7] = {8.37699/*0.848468*/, 4.02343, 92.6734, 32.9655, 0.0628411/*(1.-0.848468)*/, 211.452, 27.4743}; // L7Correction
  double *q = parsBKG_;

  if(q[6] < 0.){
    TF1 funcGamma = TF1("funcGamma", MyFunctions::gamma, 100., 550., 3);
    funcGamma.SetParameter(0, q[1]);
    funcGamma.SetParameter(1, q[2]);
    funcGamma.SetParameter(2, q[3]);
    q[6] = funcGamma.Integral(100.,550.);
  }
  if(q[7] < 0.){
    TF1 funcLandau = TF1("funcLandau", MyFunctions::landau, 100., 550., 2);
    funcLandau.SetParameter(0, q[4]);
    funcLandau.SetParameter(1, q[5]);
    q[7] = funcLandau.Integral(100.,550.);
  }

  TF1 func = TF1("func", MyFunctions::gammaland, 100., 550., 8);
  
  func.SetParameters(q);
  func.SetParameter(1, p[4]*q[1]);
    
  if(PBKGintegral_ < 0.){
    PBKGintegral_ = func.Integral(100.,550.);
  }
  double N = 1. / PBKGintegral_;

  return N*func.Eval(p[1]);
}

double IdeogramCombLikelihoodAllJets::PBKGJES1(double *x, double *p)
{
  //double q[3] = {87.1042, 5.83934, 7.02516}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  //double q[3] = {86.7801, 5.92281, 7.25256}; // L7Correction
  double *q = parsBKGJES_;

  return PBKGJES(x, p, q);
}

//double IdeogramCombLikelihoodAllJets::PBKGJES2(double *x, double *p)
//{
//  double q[3] = {85.4051, 7.02919, 11.9709};
//  return PBKGJES(x, p, q);
//}

double IdeogramCombLikelihoodAllJets::PBKGJES(double *x, double *p, double *q)
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

IdeogramCombLikelihoodAllJets::~IdeogramCombLikelihoodAllJets()
{
  delete parsCP_;
  delete parsWP_;
  delete parsUN_;
  delete parsCPJES_;
  delete parsWPJES_;
  delete parsUNJES_;
  delete parsBKG_;
  delete parsBKGJES_;
}
