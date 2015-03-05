#include "IdeogramCombLikelihoodLeptonJets.h"

#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsCP_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsWP_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsUN_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsCPJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsWPJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_parsUNJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsCP_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsWP_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsUN_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsCPJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsWPJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_parsUNJES_(0);

std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_massOffset_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_massSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_massSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_massSlopeMassJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_jesOffset_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_jesSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_jesSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::muo_jesSlopeMassJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_massOffset_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_massSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_massSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_massSlopeMassJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_jesOffset_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_jesSlopeMass_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_jesSlopeJES_(0);
std::vector<double> IdeogramCombLikelihoodLeptonJets::ele_jesSlopeMassJES_(0);

double IdeogramCombLikelihoodLeptonJets::muo_fCP_(-1.);
double IdeogramCombLikelihoodLeptonJets::muo_fWP_(-1.);
double IdeogramCombLikelihoodLeptonJets::muo_fUN_(-1.);
double IdeogramCombLikelihoodLeptonJets::ele_fCP_(-1.);
double IdeogramCombLikelihoodLeptonJets::ele_fWP_(-1.);
double IdeogramCombLikelihoodLeptonJets::ele_fUN_(-1.);

IdeogramCombLikelihoodLeptonJets::ScanPointMap IdeogramCombLikelihoodLeptonJets::PWPnormalizations_;
IdeogramCombLikelihoodLeptonJets::ScanPointMap IdeogramCombLikelihoodLeptonJets::PUNnormalizations_;

IdeogramCombLikelihoodLeptonJets::IdeogramCombLikelihoodLeptonJets()
{
  //read in electron information as well

  // parameters for mTop correct permutations
  if(!muo_parsCP_.size())
    muo_parsCP_ = readParameters("templates.muo_parsCP");
  if(!ele_parsCP_.size())
    ele_parsCP_ = readParameters("templates.ele_parsCP");
  if(ele_parsCP_.size()!=muo_parsCP_.size())std::cerr << "Error: Muon and electron config for parsCP does not match. You have to check this!" << std::endl;

  // parameters for mTop wrong parmutations
  if(!muo_parsWP_.size())
    muo_parsWP_ = readParameters("templates.muo_parsWP");
  if(!ele_parsWP_.size())
    ele_parsWP_ = readParameters("templates.ele_parsWP");
  if(ele_parsWP_.size()!=muo_parsWP_.size())std::cerr << "Error: Muon and electron config for parsWP does not match. You have to check this!" << std::endl;

  // parameters for mTop unmatched permutations
  if(!muo_parsUN_.size())
    muo_parsUN_ = readParameters("templates.muo_parsUN");
  if(!ele_parsUN_.size())
    ele_parsUN_ = readParameters("templates.ele_parsUN");
  if(ele_parsUN_.size()!=muo_parsUN_.size())std::cerr << "Error: Muon and electron config for parsUN does not match. You have to check this!" << std::endl;

  // parameters for JES correct permutations
  if(!muo_parsCPJES_.size())
    muo_parsCPJES_ = readParameters("templates.muo_parsCPJES");
  if(!ele_parsCPJES_.size())
    ele_parsCPJES_ = readParameters("templates.ele_parsCPJES");
  if(ele_parsCPJES_.size()!=muo_parsCPJES_.size())std::cerr << "Error: Muon and electron config for parsCPJES does not match. You have to check this!" << std::endl;

  // parameters for JES wrong permutations
  if(!muo_parsWPJES_.size())
    muo_parsWPJES_ = readParameters("templates.muo_parsWPJES");
  if(!ele_parsWPJES_.size())
    ele_parsWPJES_ = readParameters("templates.ele_parsWPJES");
  if(ele_parsWPJES_.size()!=muo_parsWPJES_.size())std::cerr << "Error: Muon and electron config for parsWPJES does not match. You have to check this!" << std::endl;

  // parameters for JES unmatched permutations
  if(!muo_parsUNJES_.size())
    muo_parsUNJES_ = readParameters("templates.muo_parsUNJES");
  if(!ele_parsUNJES_.size())
    ele_parsUNJES_ = readParameters("templates.ele_parsUNJES");
  if(ele_parsUNJES_.size()!=muo_parsUNJES_.size())std::cerr << "Error: Muon and electron config for parsUNJES does not match. You have to check this!" << std::endl;

  // read calibration
  if(!muo_massOffset_.size())
    muo_massOffset_ = readParameters("calibration.muo_massOffset");
  if(!muo_massSlopeMass_.size()){
    muo_massSlopeMass_ = readParameters("calibration.muo_massSlopeMass");
    if(muo_massSlopeMass_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_massSlopeMass_.size() = " << muo_massSlopeMass_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_massOffset_.size())
    ele_massOffset_ = readParameters("calibration.ele_massOffset");
  if(!ele_massSlopeMass_.size()){
    ele_massSlopeMass_ = readParameters("calibration.ele_massSlopeMass");
    if(ele_massSlopeMass_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_massSlopeMass_.size() = " << ele_massSlopeMass_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_massSlopeJES_.size()){
    muo_massSlopeJES_ = readParameters("calibration.muo_massSlopeJES");
    if(muo_massSlopeJES_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_massSlopeJES_.size() = " << muo_massSlopeJES_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_massSlopeJES_.size()){
    ele_massSlopeJES_ = readParameters("calibration.ele_massSlopeJES");
    if(ele_massSlopeJES_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_massSlopeJES_.size() = " << ele_massSlopeJES_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_massSlopeMassJES_.size()){
    muo_massSlopeMassJES_ = readParameters("calibration.muo_massSlopeMassJES");
    if(muo_massSlopeMassJES_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_massSlopeMassJES_.size() = " << muo_massSlopeMassJES_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_massSlopeMassJES_.size()){
    ele_massSlopeMassJES_ = readParameters("calibration.ele_massSlopeMassJES");
    if(ele_massSlopeMassJES_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_massSlopeMassJES_.size() = " << ele_massSlopeMassJES_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_jesOffset_.size()){
    muo_jesOffset_ = readParameters("calibration.muo_jesOffset");
    if(muo_jesOffset_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_jesOffset_.size() = " << muo_jesOffset_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_jesOffset_.size()){
    ele_jesOffset_ = readParameters("calibration.ele_jesOffset");
    if(ele_jesOffset_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_jesOffset_.size() = " << ele_jesOffset_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_jesSlopeMass_.size()){
    muo_jesSlopeMass_ = readParameters("calibration.muo_jesSlopeMass");
    if(muo_jesSlopeMass_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_jesSlopeMass_.size() = " << muo_jesSlopeMass_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_jesSlopeMass_.size()){
    ele_jesSlopeMass_ = readParameters("calibration.ele_jesSlopeMass");
    if(ele_jesSlopeMass_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_jesSlopeMass_.size() = " << ele_jesSlopeMass_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_jesSlopeJES_.size()){
    muo_jesSlopeJES_ = readParameters("calibration.muo_jesSlopeJES");
    if(muo_jesSlopeJES_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_jesSlopeJES_.size() = " << muo_jesSlopeJES_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_jesSlopeJES_.size()){
    ele_jesSlopeJES_ = readParameters("calibration.ele_jesSlopeJES");
    if(ele_jesSlopeJES_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_jesSlopeJES_.size() = " << ele_jesSlopeJES_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(!muo_jesSlopeMassJES_.size()){
    muo_jesSlopeMassJES_ = readParameters("calibration.muo_jesSlopeMassJES");
    if(muo_jesSlopeMassJES_.size() != muo_massOffset_.size()){
      std::cout << "Different number of calibration constants! muo_jesSlopeMassJES_.size() = " << muo_jesSlopeMassJES_.size() << ", muo_massOffset_.size() = " << muo_massOffset_.size() << std::endl;
    }
  }
  if(!ele_jesSlopeMassJES_.size()){
    ele_jesSlopeMassJES_ = readParameters("calibration.ele_jesSlopeMassJES");
    if(ele_jesSlopeMassJES_.size() != ele_massOffset_.size()){
      std::cout << "Different number of calibration constants! ele_jesSlopeMassJES_.size() = " << ele_jesSlopeMassJES_.size() << ", ele_massOffset_.size() = " << ele_massOffset_.size() << std::endl;
    }
  }

  if(muo_fCP_  == -1) muo_fCP_  = po::GetOption<double>("templates.muo_fCP");
  if(muo_fWP_  == -1) muo_fWP_  = po::GetOption<double>("templates.muo_fWP");
  if(muo_fUN_  == -1) muo_fUN_  = po::GetOption<double>("templates.muo_fUN");
  if(ele_fCP_  == -1) ele_fCP_  = po::GetOption<double>("templates.ele_fCP");
  if(ele_fWP_  == -1) ele_fWP_  = po::GetOption<double>("templates.ele_fWP");
  if(ele_fUN_  == -1) ele_fUN_  = po::GetOption<double>("templates.ele_fUN");


  if(muo_parsCP_   .size()!=12||
     muo_parsWP_   .size()!=12||
     muo_parsUN_   .size()!=12||
     muo_parsCPJES_.size()!=12||
     muo_parsWPJES_.size()!=12||
     muo_parsUNJES_.size()!=12||
     ele_parsCP_   .size()!=12||
     ele_parsWP_   .size()!=12||
     ele_parsUN_   .size()!=12||
     ele_parsCPJES_.size()!=12||
     ele_parsWPJES_.size()!=12||
     ele_parsUNJES_.size()!=12)std::cout << "Error: At least one of the parameter configs for the templates has too many or too few parameters (checked for 12). Check, please!" << std::endl;


};

double IdeogramCombLikelihoodLeptonJets::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool useCalib = true;
  
  //* muons from config
  double fCP = muo_fCP_; // combination type 1
  double fWP = muo_fWP_; // combination types 2 3 4
  double fUN = muo_fUN_; // combination types -2 -1 6
  //*/
  
  // for IdeogramMinimizer
  if (useFixedParams_) {
    double ep[9];
    for (int i = 0; i < 9; ++i) {
      ep[i] = fp_[i];
    }
    p = ep;
  }
  
  //* electrons from config
  if (p[3] == 11) {
    fCP = ele_fCP_;
    fWP = ele_fWP_;
    fUN = ele_fUN_;
  }
  //*
  
  /* Fall11 P11
     double fCP = 0.415661;
     double fWP = 0.201059;
     double fUN = 0.38328;
  //*/
  
  /* worst case, improvable by f(w_i)
     double fCP = 0.439252;
     double fWP = 0.211647;
     double fUN = 0.349102;
  //*/
  
  if (onlyCP) {
    fWP = 0;
    fUN = 0;
  }
  
  
  double MassOffset       = muo_massOffset_       .at(0);
  double MassSlopeMass    = muo_massSlopeMass_    .at(0);
  double MassSlopeJES     = muo_massSlopeJES_     .at(0);
  double MassSlopeMassJES = muo_massSlopeMassJES_ .at(0);

  double JESOffset        = muo_jesOffset_        .at(0);
  double JESSlopeMass     = muo_jesSlopeMass_     .at(0);
  double JESSlopeJES      = muo_jesSlopeJES_      .at(0);
  double JESSlopeMassJES  = muo_jesSlopeMassJES_  .at(0);
  
  double fastMassOffset       = muo_massOffset_       .at(1);
  double fastMassSlopeMass    = muo_massSlopeMass_    .at(1);
  double fastMassSlopeJES     = muo_massSlopeJES_     .at(1);
  double fastMassSlopeMassJES = muo_massSlopeMassJES_ .at(1);

  double fastJESOffset        = muo_jesOffset_        .at(1);
  double fastJESSlopeMass     = muo_jesSlopeMass_     .at(1);
  double fastJESSlopeJES      = muo_jesSlopeJES_      .at(1);
  double fastJESSlopeMassJES  = muo_jesSlopeMassJES_  .at(1);

  if (p[3] == 11) {
    MassOffset       = ele_massOffset_       .at(0);
    MassSlopeMass    = ele_massSlopeMass_    .at(0);
    MassSlopeJES     = ele_massSlopeJES_     .at(0);
    MassSlopeMassJES = ele_massSlopeMassJES_ .at(0);

    JESOffset        = ele_jesOffset_        .at(0);
    JESSlopeMass     = ele_jesSlopeMass_     .at(0);
    JESSlopeJES      = ele_jesSlopeJES_      .at(0);
    JESSlopeMassJES  = ele_jesSlopeMassJES_  .at(0);
    
    fastMassOffset       = ele_massOffset_       .at(1);
    fastMassSlopeMass    = ele_massSlopeMass_    .at(1);
    fastMassSlopeJES     = ele_massSlopeJES_     .at(1);
    fastMassSlopeMassJES = ele_massSlopeMassJES_ .at(1);

    fastJESOffset        = ele_jesOffset_        .at(1);
    fastJESSlopeMass     = ele_jesSlopeMass_     .at(1);
    fastJESSlopeJES      = ele_jesSlopeJES_      .at(1);
    fastJESSlopeMassJES  = ele_jesSlopeMassJES_  .at(1);
  }
  
  double m = x[0];
  double j = x[1];
  
  if (useCalib) {
    x[1] = j + JESOffset + JESSlopeMass*(m-172.5) + JESSlopeJES*(j-1.) + JESSlopeMassJES*(m-172.5)*(j-1.);
    x[0] = m + MassOffset + MassSlopeMass*(m-172.5) + MassSlopeJES*(j-1.) + MassSlopeMassJES*(m-172.5)*(j-1.);
  }
  
  if (p[7] == 1) {
    x[1] = j + fastJESOffset + fastJESSlopeMass*(m-172.5) + fastJESSlopeJES*(j-1.) + fastJESSlopeMassJES*(m-172.5)*(j-1.);
    x[0] = m + fastMassOffset + fastMassSlopeMass*(m-172.5) + fastMassSlopeJES*(j-1.) + fastMassSlopeMassJES*(m-172.5)*(j-1.);
  }
  
  //return p[0] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[0] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));
}


double IdeogramCombLikelihoodLeptonJets::PCP(double *x, double *p) {
  double* q = &muo_parsCP_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsCP_[0];

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


namespace cb {
  double A(double alpha, double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(double alpha, double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}


double IdeogramCombLikelihoodLeptonJets::PWP(double* x, double* p)
{
  double* q = &muo_parsWP_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsWP_[0];
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double power  =  15;
  double t = (p[1] - mu) / sigma;
  
  ScanPoint currentPoint = std::make_tuple(p[3],x[0],x[1]);
  ScanPointMap::const_iterator iter = PWPnormalizations_.find(currentPoint);
  double currentNormalization = -1.;
  if(iter == PWPnormalizations_.end()){
    double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
    double N2 = cb::A(alpha,power) / (-1.+power) * (
						    (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
						    -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
						    );
    currentNormalization = N1 + N2;
    PWPnormalizations_[currentPoint] = currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  
  if(t < alpha)
    return 1./currentNormalization * TMath::Exp(-t*t/2);
  else
    return 1./currentNormalization * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}


double IdeogramCombLikelihoodLeptonJets::PUN(double* x, double* p)
{
  double* q = &muo_parsUN_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsUN_[0];
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
  double power  =  5;
  double t = (p[1] - mu) / sigma;
  
  ScanPoint currentPoint = std::make_tuple(p[3],x[0],x[1]);
  ScanPointMap::const_iterator iter = PUNnormalizations_.find(currentPoint);
  double currentNormalization = -1.;
  if(iter == PUNnormalizations_.end()){
    double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
    double N2 = cb::A(alpha,power) / (-1.+power) * (
						    (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
						    -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
						    );
    currentNormalization = N1 + N2;
    PUNnormalizations_[currentPoint] = currentNormalization;
  }
  else{
    currentNormalization = iter->second;
  }
  
  if(t < alpha)
    return 1./currentNormalization * TMath::Exp(-t*t/2);
  else
    return 1./currentNormalization * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}

double IdeogramCombLikelihoodLeptonJets::PCPJES(double* x, double* p)
{
  double* q = &muo_parsCPJES_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsCPJES_[0];
  
  //* W Mass
  //std::cout << p[4] << " ";
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
   
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
     double mu       = 1.36597e-02 - 1.00417e-01 * (x[1]-1.);
     double sigma    = 4.41279e-02 + 3.22201e-02 * (x[1]-1.);
     double width    = 6.67793e-02 - 4.12224e-02 * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
	
  /* Constraint
     return TMath::Gaus(x[1], 1, 0.1, kTRUE);
  //*/
	
  /* b scale estimator
     double mu       = 1.23115e+00 - 4.21823e-03 * (x[0]-172.5) + (-1.48697e+00 - 4.92782e-03 * (x[0]-172.5)) * (x[1]-1.);
     double sigma    = 0.001;
     double width    = 1.29952e+00 + 1.10928e-02 * (x[0]-172.5) + (-1.14917e+00 - 3.74889e-02 * (x[0]-172.5)) * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}

double IdeogramCombLikelihoodLeptonJets::PWPJES(double* x, double* p)
{
  double* q = &muo_parsWPJES_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsWPJES_[0];
  
  //* W Mass
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
      
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
     double mu       = 1.94840e-02 - 5.30632e-02 * (x[1]-1.);
     double sigma    = 5.86415e-02 + 2.58850e-02 * (x[1]-1.);
     double width    = 1.22249e-01 - 1.36272e-01 * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
  
  /* b scale estimator
     double mu       = 7.18283e-02 + 1.68825e-01 * (x[1]-1.);
     double sigma    = 9.92390e-01 - 8.84443e-01 * (x[1]-1.);
     double width    = 2.15202e+00 - 2.33270e+00 * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}

double IdeogramCombLikelihoodLeptonJets::PUNJES(double* x, double* p)
{
  double* q = &muo_parsUNJES_[0]; //could use it as vector (range check)
  if (p[3] == 11) q = &ele_parsUNJES_[0];

  //* W Mass
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
    
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
     double mu       = 9.99882e-03 + 2.33525e-02 * (x[1]-1.);
     double sigma    = 5.51757e-02 - 2.62315e-02 * (x[1]-1.);
     double width    = 1.35272e-01 - 3.44117e-03 * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
  
  /* b scale estimator
     double mu       = 9.45273e-01 - 1.07859e+00 * (x[1]-1.) - 0.2;
     double sigma    = 0.001;
     double width    = 2.14746e+00 - 2.29404e+00 * (x[1]-1.);
  
     return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}
