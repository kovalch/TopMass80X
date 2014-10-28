#include "IdeogramCombLikelihoodAllJets.h"

#include "ProgramOptionsReader.h"

#include <stdlib.h>

#include "TF1.h"

typedef ProgramOptionsReader po;

std::vector<double> IdeogramCombLikelihoodAllJets::parsCP_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsWP_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsUN_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsCPJES_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsWPJES_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsUNJES_(0);

std::vector<double> IdeogramCombLikelihoodAllJets::massOffset_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::massSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::massSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::massSlopeMassJES_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::jesOffset_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::jesSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::jesSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::jesSlopeMassJES_(0);

std::vector<double> IdeogramCombLikelihoodAllJets::parsBKG_(0);
std::vector<double> IdeogramCombLikelihoodAllJets::parsBKGJES_(0);

double IdeogramCombLikelihoodAllJets::fCP_(-1);
double IdeogramCombLikelihoodAllJets::fWP_(-1);
double IdeogramCombLikelihoodAllJets::fUN_(-1);

double IdeogramCombLikelihoodAllJets::IntegrationRangeMax_(215.);
double IdeogramCombLikelihoodAllJets::PBKGintegral_(-1);
IdeogramCombLikelihoodAllJets::ScanPointMap IdeogramCombLikelihoodAllJets::PWPnormalizations_;
IdeogramCombLikelihoodAllJets::ScanPointMap IdeogramCombLikelihoodAllJets::PUNnormalizations_;

IdeogramCombLikelihoodAllJets::IdeogramCombLikelihoodAllJets() :
  fSig_(-1.)
  {
  bool onlyCP = false;

  // parameters for mTop correct permutations
  if(!parsCP_.size())
    parsCP_ = readParameters("templates.parsCP");

  // parameters for mTop wrong parmutations
  if(!parsWP_.size())
    parsWP_ = readParameters("templates.parsWP");

  // parameters for mTop unmatched permutations
  if(!parsUN_.size())
    parsUN_ = readParameters("templates.parsUN");

  // parameters for JES correct permutations
  if(!parsCPJES_.size())
    parsCPJES_ = readParameters("templates.parsCPJES");

  // parameters for JES wrong permutations
  if(!parsWPJES_.size())
    parsWPJES_ = readParameters("templates.parsWPJES");

  // parameters for JES unmatched permutations
  if(!parsUNJES_.size())
    parsUNJES_ = readParameters("templates.parsUNJES");

  // read calibration
  if(!massOffset_.size())
    massOffset_ = readParameters("calibration.massOffset");
  if(!massSlopeMass_.size()){
    massSlopeMass_ = readParameters("calibration.massSlopeMass");
    if(massSlopeMass_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! massSlopeMass_.size() = " << massSlopeMass_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!massSlopeJES_.size()){
    massSlopeJES_ = readParameters("calibration.massSlopeJES");
    if(massSlopeJES_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! massSlopeJES_.size() = " << massSlopeJES_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!massSlopeMassJES_.size()){
    massSlopeMassJES_ = readParameters("calibration.massSlopeMassJES");
    if(massSlopeMassJES_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! massSlopeMassJES_.size() = " << massSlopeMassJES_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!jesOffset_.size()){
    jesOffset_ = readParameters("calibration.jesOffset");
    if(jesOffset_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! jesOffset_.size() = " << jesOffset_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!jesSlopeMass_.size()){
    jesSlopeMass_ = readParameters("calibration.jesSlopeMass");
    if(jesSlopeMass_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! jesSlopeMass_.size() = " << jesSlopeMass_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!jesSlopeJES_.size()){
    jesSlopeJES_ = readParameters("calibration.jesSlopeJES");
    if(jesSlopeJES_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! jesSlopeJES_.size() = " << jesSlopeJES_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }
  if(!jesSlopeMassJES_.size()){
    jesSlopeMassJES_ = readParameters("calibration.jesSlopeMassJES");
    if(jesSlopeMassJES_.size() != massOffset_.size()){
      std::cout << "Different number of calibration constants! jesSlopeMassJES_.size() = " << jesSlopeMassJES_.size() << ", massOffset_.size() = " << massOffset_.size() << std::endl;
    }
  }

  // parameters for mTop background
  if(!parsBKG_.size()){
    parsBKG_ = readParameters("templates.parsBKG");
    parsBKG_.push_back(-1);
    parsBKG_.push_back(-1);
  }

  // parameters for JES background
  if(!parsBKGJES_.size())
    parsBKGJES_ = readParameters("templates.parsBKGJES");
 
  if(fCP_ == -1) fCP_ = po::GetOption<double>("templates.fCP");
  if(fWP_ == -1) fWP_ = po::GetOption<double>("templates.fWP");
  if(fUN_ == -1) fUN_ = po::GetOption<double>("templates.fUN");
  
  if (onlyCP) {
    fCP_ = 1;
    fWP_ = 0;
    fUN_ = 0;
  }
}

double IdeogramCombLikelihoodAllJets::Evaluate(double *x, double *p) {
  bool useCalib   = true;
  bool onlySIG    = false;
  //bool fSigOfMtop = false;
  
  // for IdeogramMinimizer
  if (useFixedParams_) {
    double ep[9];
    for (int i = 0; i < 9; ++i) {
      ep[i] = fp_[i];
    }
    p = ep;
  }
  
  if(p[6]){
    fCP_ -= 0.5*p[6];
    fWP_ -= 0.5*p[6];
    fUN_ += p[6];
  }

  //if(fSigOfMtop){
  //  double a0    = -0.00464347142282148;   // CSVT, P(chi2) weighted
  //  double a1    =  0.0000348026883089243; // CSVT, P(chi2) weighted
  //  double b0    =  477200000000; // <- Ahrens //  497700000000; // <- Kidonakis //  501000000000; // <- Langenfeld //
  //  double b1    = -3322000000;  // <- Ahrens //  -3375000000;  // <- Kidonakis //  -3385000000;  // <- Langenfeld //
  //  double b2    =  10050000;   // <- Ahrens //    9625000;    // <- Kidonakis //    9553000;    // <- Langenfeld //
  //  double b3    = -12340;     // <- Ahrens //    -10770;     // <- Kidonakis //    -10470;     // <- Langenfeld //
  //  double Lint  = 3544.844;
  //  double Ndata = 2418.*0.437638;
  //
  //  fSig_ = (a0+a1*x[0])*(b0+b1*x[0]+b2*pow(x[0],2)+b3*pow(x[0],3))/pow(x[0],4)*Lint/Ndata;
  //}

  if (useCalib) {

    for(unsigned int i = 0, l = massOffset_.size(); i < l; ++i){
      // skip last calibration step if not using FastSim
      if(i == l-1 && p[7] == 0) break;
      double m0 = x[0];
      double j0 = x[1];

      x[0] = m0 + massOffset_[i] + massSlopeMass_[i]*(m0-172.5) + massSlopeJES_[i]*(j0-1.) + massSlopeMassJES_[i]*(m0-172.5)*(j0-1.);
      x[1] = j0 +  jesOffset_[i] +  jesSlopeMass_[i]*(m0-172.5) +  jesSlopeJES_[i]*(j0-1.) +  jesSlopeMassJES_[i]*(m0-172.5)*(j0-1.);
    }
  }

  if(onlySIG){
    fSig_ = 1.;
  }
  else{
    fSig_ = x[2];
  }

  if(fSig_ == -1){
    fSig_ = po::GetOption<double>("templates.fSig");
  }

  fCP_ = x[3];
  if(!fUN_)
    fWP_ = 1.0-x[3];
  else{
    static double rfWP = fWP_/(fUN_+fWP_);
    static double rfUN = fUN_/(fUN_+fWP_);
    fWP_ = rfWP*(1.0-x[3]);
    fUN_ = rfUN*(1.0-x[3]);
  }

  return p[0] * (fSig_ * (fCP_ * PCP (x,p) * (PCPJES1 (x,p))
		       +  fWP_ * PWP (x,p) * (PWPJES1 (x,p))
     +  ((fUN_==0) ? 0 : (fUN_ * PUN (x,p) * (PUNJES1 (x,p)))) )
	       + (1.-fSig_)    * PBKG(x,p) * (PBKGJES1(x,p)));
}


double IdeogramCombLikelihoodAllJets::PCP(double *x, double *p) {
  const std::vector<double> &q = parsCP_;

  double mu    = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
 }

double IdeogramCombLikelihoodAllJets::PWP(double* x, double* p)
{
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
    //TF1 func = TF1("func", langau, 100., IntegrationRangeMax_, 5);
    TF1 func = TF1("func", lanvog, 100., IntegrationRangeMax_, 6);
    func.SetParameters(funcPars);
    currentNormalization = 1. / func.Integral(100.,IntegrationRangeMax_);
    PWPnormalizations_[currentPoint] = currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  return currentNormalization*lanvog(&p[1],&funcPars[0]);
}

double IdeogramCombLikelihoodAllJets::PUN(double* x, double* p)
{
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
    //TF1 func = TF1("func", langau, 100., IntegrationRangeMax_, 5);
    TF1 func = TF1("func", lanvog, 100., IntegrationRangeMax_, 6);
    func.SetParameters(funcPars);
    currentNormalization = 1. / func.Integral(100.,IntegrationRangeMax_);
    PUNnormalizations_[currentPoint] = currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  return currentNormalization*lanvog(&p[1],&funcPars[0]);
}

double IdeogramCombLikelihoodAllJets::PCPJES1(double *x, double *p)
{
  const std::vector<double> &q = parsCPJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PWPJES1(double *x, double *p)
{
  const std::vector<double> &q = parsWPJES_;

  return PJESOTHER(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PUNJES1(double *x, double *p)
{
  const std::vector<double> &q = parsUNJES_;

  return PJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PJES(double* x, double* p, const std::vector<double> &q)
{
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2]  + q[3]  * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6]  + q[7]  * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
   
  double t1 = (p[2] - mu) / sigma1;
  double t2 = (p[2] - mu) / sigma2;
  
  static const double sqrt2  = sqrt(2.0);
  double N = sqrt2/(sigma1+sigma2);
  
  if(t1 < 0)
    return N * TMath::Exp(-t1*t1/2);
  else
    return N * TMath::Exp(-t2*t2/2);
}

double IdeogramCombLikelihoodAllJets::PJESOTHER(double* x, double* p, const std::vector<double> &q)
{
  double mu     = q[ 0] + q[ 1] * (x[0]-172.5) + (q[ 2] + q[ 3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[ 4] + q[ 5] * (x[0]-172.5) + (q[ 6] + q[ 7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[ 8] + q[ 9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double sigma3 = q[12] + q[13] * (x[0]-172.5) + (q[14] + q[15] * (x[0]-172.5)) * (x[1]-1.);
  
  double shift1 = 7.5;
  double shift3 = 7.5;
  double norm1 = 0.25;
  double norm2 = 0.5;
  if(q.size() == 22){
    shift1 = q[16];
    shift3 = q[17];
    norm1  = q[18];
    norm2  = q[19];
  }
  double norm3 = 1-norm1-norm2;

  double t1 = (p[2] - (mu - shift1)) / sigma1;
  double t2 = (p[2] -  mu          ) / sigma2;
  double t3 = (p[2] - (mu + shift3)) / sigma3;
  
  static const double sqrt2 = sqrt(2.0);
  double N1 = norm1/(sqrt2*sigma1);
  double N2 = norm2/(sqrt2*sigma2);
  double N3 = norm3/(sqrt2*sigma3);

  return ((N1 * TMath::Exp(-t1*t1/2)) + (N2 * TMath::Exp(-t2*t2/2)) + (N3 * TMath::Exp(-t3*t3/2)));
}

double IdeogramCombLikelihoodAllJets::PBKG(double *x, double *p)
{
  // [0]*TMath::GammaDist(x,[1],[2],[3])+[4]*TMath::Landau(x,[5],[6])
  std::vector<double> &q = parsBKG_;

  if(q[6] < 0.){
    TF1 funcGamma = TF1("funcGamma", gamma, 100., IntegrationRangeMax_, 3);
    funcGamma.SetParameter(0, p[4]>0 ? p[4]*q[1] : q[1]);
    funcGamma.SetParameter(1, q[2]);
    funcGamma.SetParameter(2, q[3]);
    q[6] = funcGamma.Integral(100.,IntegrationRangeMax_);
  }
  if(q[7] < 0.){
    TF1 funcLandau = TF1("funcLandau", landau, 100., IntegrationRangeMax_, 2);
    funcLandau.SetParameter(0, q[4]);
    funcLandau.SetParameter(1, q[5]);
    q[7] = funcLandau.Integral(100.,IntegrationRangeMax_);
  }

  if(PBKGintegral_ < 0.){
    TF1 func = TF1("func", gammaland, 100., IntegrationRangeMax_, 8);
    q[1] = p[4]>0 ? p[4]*q[1] : q[1];
    func.SetParameters(&q[0]);
    PBKGintegral_ = 1. / func.Integral(100.,IntegrationRangeMax_);
  }

  return PBKGintegral_*gammaland(&p[1],&q[0]);
}

double IdeogramCombLikelihoodAllJets::PBKGJES1(double *x, double *p)
{
  const std::vector<double> &q = parsBKGJES_;

  return PBKGJES(x, p, q);
}

double IdeogramCombLikelihoodAllJets::PBKGJES(double *x, double *p, const std::vector<double> &q)
{
  double mu     = p[5]>0 ? q[0]*p[5] : q[0];
  double sigma1 = q[1];
  double sigma2 = q[2];
   
  double t1 = (p[2] - mu) / sigma1;
  double t2 = (p[2] - mu) / sigma2;
  
  static const double sqrt2 = sqrt(2);
  double N = sqrt2/(sigma1+sigma2);
  
  if(t1 < 0)
    return N * TMath::Exp(-t1*t1/2);
  else
    return N * TMath::Exp(-t2*t2/2);
}

IdeogramCombLikelihoodAllJets::~IdeogramCombLikelihoodAllJets()
{
}
