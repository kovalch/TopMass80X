#include "IdeogramCombLikelihoodAllJets.h"

#include "ProgramOptionsReader.h"

#include <stdlib.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "TF1.h"

typedef ProgramOptionsReader po;

IdeogramCombLikelihoodAllJets::IdeogramCombLikelihoodAllJets() :
  parsBKG_  (0), parsBKGJES_(0),
  fSig_(-1.),
  PBKGintegral_(-1)
  {

  // parameters for mTop background
  if(!parsBKG_.size()){
    parsBKG_ = readParameters("templates.parsBKG");
    parsBKG_.push_back(-1);
    parsBKG_.push_back(-1);
  }

  // parameters for JES background
  if(!parsBKGJES_.size())
    parsBKGJES_ = readParameters("templates.parsBKGJES");

  if(fSig_ == -1) fSig_ = po::GetOption<double>("templates.fSig");

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

    for(unsigned int i = 0, l = massOffset_.size(); i < l; ++i){
      // skip last calibration step if not using FastSim
      if(i == l-1 && p[6] == 0) break;
      double m0 = x[0];
      double j0 = x[1];

      x[0] = m0 + massOffset_[i] + massSlopeMass_[i]*(m0-172.5) + massSlopeJES_[i]*(j0-1.) + massSlopeMassJES_[i]*(m0-172.5)*(j0-1.);
      x[1] = j0 +  jesOffset_[i] +  jesSlopeMass_[i]*(m0-172.5) +  jesSlopeJES_[i]*(j0-1.) +  jesSlopeMassJES_[i]*(m0-172.5)*(j0-1.);
    }

    //double MassOffset0       = 0.; // 3.31715e-01;
    //double MassSlopeMass0    = 0.; // 8.22606e-03;
    //double MassSlopeJES0     = 0.; // 4.60359e-01;
    //double MassSlopeMassJES0 = 0.; //-2.22772e-01;
    //
    //double JESOffset0        = 0.; //-8.78984e-04;
    //double JESSlopeMass0     = 0.; //-4.66983e-05;
    //double JESSlopeJES0      = 0.; //-5.75365e-02;
    //double JESSlopeMassJES0  = 0.; // 1.56276e-03;
    //
    //double m0 = x[0];
    //double j0 = x[1];
    //
    //double j1 = j0 +  JESOffset0 +  JESSlopeMass0*(m0-172.5) +  JESSlopeJES0*(j0-1.) +  JESSlopeMassJES0*(m0-172.5)*(j0-1.);
    //double m1 = m0 + MassOffset0 + MassSlopeMass0*(m0-172.5) + MassSlopeJES0*(j0-1.) + MassSlopeMassJES0*(m0-172.5)*(j0-1.);
    //
    //// do additional FastSim calibration if p[6] == 1 (command line parameter -F 1)
    //double MassOffset1       = (p[6] == 1) ?  1.83190 : 0.;
    //double MassSlopeMass1    = 0.;
    //double MassSlopeJES1     = (p[6] == 1) ? -5.36113 : 0.;
    //double MassSlopeMassJES1 = 0.;
	//
    //double JESOffset1        = (p[6] == 1) ? -2.11957e-02 : 0.;
    //double JESSlopeMass1     = 0.;
    //double JESSlopeJES1      = (p[6] == 1) ?  3.15759e-02 : 0.;
    //double JESSlopeMassJES1  = 0.;
    //
    //x[1] = j1 +  JESOffset1 +  JESSlopeMass1*(m1-172.5) +  JESSlopeJES1*(j1-1.) +  JESSlopeMassJES1*(m1-172.5)*(j1-1.);
    //x[0] = m1 + MassOffset1 + MassSlopeMass1*(m1-172.5) + MassSlopeJES1*(j1-1.) + MassSlopeMassJES1*(m1-172.5)*(j1-1.);
  }
  
  return p[0] * (fSig_ * (fCP_ * PCP (x,p) * (PCPJES1 (x,p))
		               +  fWP_ * PWP (x,p) * (PWPJES1 (x,p))
		               +  fUN_ * PUN (x,p) * (PUNJES1 (x,p)))
	       + (1.-fSig_)        * PBKG(x,p) * (PBKGJES1(x,p)));
}


double IdeogramCombLikelihoodAllJets::PCP(double *x, double *p) {
  //double q[8] = {172.324, 0.986106, 85.9343, 0.810183, 7.18309, 0.0675277, 6.85536, 0.256293}; // used for measurement
  //double q[8] = {171.321, 0.96406, 82.4588, 0.667791, 7.17632, 0.0636133, 5.65153, 0.0753321}; // L7Correction
  const std::vector<double> &q = parsCP_;

  double mu    = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
 }

double IdeogramCombLikelihoodAllJets::PWP(double* x, double* p)
{
  //double q[12] = {175.941, 0.661709, 73.1694, 0.36766, 21.61, 0.100468, 6.74084, 0.561399, 7.44604, -0.0502921, -1.84577, 0.100049}; // used for measurement
  //double q[16] = {166.62, 0.479989, 65.8855, 0.817456, 172.708, 1.0786, 92.6114, -1.10235, 15.9693, -0.0223771, 0.348473, 0.177825, 7.08448, 0.115902, 6.08694, -1.16462}; // L7Correction
  const std::vector<double> &q = parsWP_;

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = q[16];

  double funcPars[] = {1., muL, muG, sigmaL, sigmaG, ratio};
    
  ScanPoint currentPoint = std::make_pair(x[0],x[1]);
  ScanPointMap::const_iterator iter = PWPnormalizations_.find(currentPoint);
  double currentNormalization = -1.;
  if(iter == PWPnormalizations_.end()){
    //TF1 func = TF1("func", langau, 100., 550., 5);
    TF1 func = TF1("func", lanvog, 100., 550., 6);
    func.SetParameters(funcPars);
    currentNormalization = func.Integral(100.,550.);
    PWPnormalizations_[currentPoint] = 1. / currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  return currentNormalization*lanvog(&p[1],&funcPars[0]);
}

double IdeogramCombLikelihoodAllJets::PUN(double* x, double* p)
{
  //double q[12] = {175.453, 0.736703, 70.2111, 1.69023, 21.1899, 0.189805, 15.3274, 1.60513, 8.87531, -0.0238031, -0.117816, 0.760937}; // used for measurement
  //double q[16] = {178.286, 0.26325, 45.3131, 0.692725, 173.347, 0.951113, 79.5679, 1.57497, 22.5683, 0.0117107, 7.02635, 0.604059, 9.93498, 0.0310997, -12.2806, 0.000101749}; // L7Correction
  const std::vector<double> &q = parsUN_;

  double muL    = q[0]  + q[1]  * (x[0]-172.5) + (q[2]   + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double muG    = q[4]  + q[5]  * (x[0]-172.5) + (q[6]   + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigmaL = q[8]  + q[9]  * (x[0]-172.5) + (q[10]  + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigmaG = q[12] + q[13] * (x[0]-172.5) + (q[14]  + q[15] * (x[0]-172.5)) * (x[1]-1.);

  double ratio  = q[16];

  double funcPars[] = {1., muL, muG, sigmaL, sigmaG, ratio};
    
  ScanPoint currentPoint = std::make_pair(x[0],x[1]);
  ScanPointMap::const_iterator iter = PUNnormalizations_.find(currentPoint);
  double currentNormalization = -1.;
  if(iter == PUNnormalizations_.end()){
    //TF1 func = TF1("func", langau, 100., 550., 5);
    TF1 func = TF1("func", lanvog, 100., 550., 6);
    func.SetParameters(funcPars);
    currentNormalization = func.Integral(100.,550.);
    PUNnormalizations_[currentPoint] = 1. / currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  return currentNormalization*lanvog(&p[1],&funcPars[0]);
}

double IdeogramCombLikelihoodAllJets::PCPJES1(double *x, double *p)
{
  //double q[12] = {84.4587, -0.00524409, 57.1528, 0.178028, 4.25284, 0.000167644, 11.7426, 0.0582509, 5.43288, -0.00549247, -6.9163, -0.0576681}; // used for measurement
  //double q[12] = {80.7172, -0.0192033, 58.9908, 0.216647, 3.56043, -0.00835824, 11.4794, 0.0100302, 5.81594, 0.00560078, -8.77514, -0.240875}; // L7Correction
  const std::vector<double> &q = parsCPJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PWPJES1(double *x, double *p)
{
  //double q[12] = {84.042, 0.0111381, 27.5838, -0.179698, 4.28814, -0.00702476, 5.27829, -0.234302, 6.28184, -0.0166415, -6.23916, 0.181574}; // used for measurement
  //double q[12] = {82.6125, 0.0113286, 14.1711, 0.94118, 4.20692, 0.000889332, -1.78333, 0.300516, 6.74883, -0.0193121, -6.19779, -0.343831}; // L7Correction
  const std::vector<double> &q = parsWPJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PUNJES1(double *x, double *p)
{
  //double q[12] = {84.6224, 0.0203645, 16.0195, -0.293667, 4.77153, 0.0196543, 2.16641, -0.150288, 7.09421, -0.00659033, -6.2366, 0.235723}; // used for measurement
  //double q[12] = {83.069, 0.0126392, 15.8747, 0.101566, 4.56234, 0.00662185, 0.145635, 0.188651, 7.54043, -0.00766768, -10.0228, 0.00347318}; // L7Correction
  const std::vector<double> &q = parsUNJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PJES(double* x, double* p, const std::vector<double> &q)
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
  std::vector<double> &q = parsBKG_;

  if(q[6] < 0.){
    TF1 funcGamma = TF1("funcGamma", gamma, 100., 550., 3);
    funcGamma.SetParameter(0, q[1]);
    funcGamma.SetParameter(1, q[2]);
    funcGamma.SetParameter(2, q[3]);
    q[6] = funcGamma.Integral(100.,550.);
  }
  if(q[7] < 0.){
    TF1 funcLandau = TF1("funcLandau", landau, 100., 550., 2);
    funcLandau.SetParameter(0, q[4]);
    funcLandau.SetParameter(1, q[5]);
    q[7] = funcLandau.Integral(100.,550.);
  }

  if(PBKGintegral_ < 0.){
    TF1 func = TF1("func", gammaland, 100., 550., 8);
    q.at(1) = p[4]*q.at(1);
    func.SetParameters(&q.at(0));
    PBKGintegral_ = 1. / func.Integral(100.,550.);
  }

  return PBKGintegral_*gammaland(&p[1],&q.at(0));
}

double IdeogramCombLikelihoodAllJets::PBKGJES1(double *x, double *p)
{
  //double q[3] = {87.1042, 5.83934, 7.02516}; // prob > 0.09 & dRbb > 1.5 - mean top mass, CSVT
  //double q[3] = {86.7801, 5.92281, 7.25256}; // L7Correction
  const std::vector<double> &q = parsBKGJES_;

  return PBKGJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PBKGJES(double *x, double *p, const std::vector<double> &q)
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
}
