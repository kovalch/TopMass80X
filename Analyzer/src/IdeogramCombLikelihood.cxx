#include "IdeogramCombLikelihood.h"

#include "ProgramOptionsReader.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

std::vector<double> IdeogramCombLikelihood::parsCP_(0);
std::vector<double> IdeogramCombLikelihood::parsWP_(0);
std::vector<double> IdeogramCombLikelihood::parsUN_(0);
std::vector<double> IdeogramCombLikelihood::parsCPJES_(0);
std::vector<double> IdeogramCombLikelihood::parsWPJES_(0);
std::vector<double> IdeogramCombLikelihood::parsUNJES_(0);

std::vector<double> IdeogramCombLikelihood::massOffset_(0);
std::vector<double> IdeogramCombLikelihood::massSlopeMass_(0);
std::vector<double> IdeogramCombLikelihood::massSlopeJES_(0);
std::vector<double> IdeogramCombLikelihood::massSlopeMassJES_(0);
std::vector<double> IdeogramCombLikelihood::jesOffset_(0);
std::vector<double> IdeogramCombLikelihood::jesSlopeMass_(0);
std::vector<double> IdeogramCombLikelihood::jesSlopeJES_(0);
std::vector<double> IdeogramCombLikelihood::jesSlopeMassJES_(0);

double IdeogramCombLikelihood::fCP_(-1);
double IdeogramCombLikelihood::fWP_(-1);
double IdeogramCombLikelihood::fUN_(-1);

IdeogramCombLikelihood::IdeogramCombLikelihood():
		//parsCP_   (0), parsWP_    (0), parsUN_   (0),
		//parsCPJES_(0), parsWPJES_ (0), parsUNJES_(0),
		//massOffset_(0), massSlopeMass_(0), massSlopeJES_(0), massSlopeMassJES_(0),
		//jesOffset_ (0), jesSlopeMass_ (0), jesSlopeJES_ (0), jesSlopeMassJES_ (0),
		//fCP_(-1.), fWP_(-1.), fUN_(-1.),
		useFixedParams_(false)
{
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

  if(fCP_ == -1) fCP_ = po::GetOption<double>("templates.fCP");
  if(fWP_ == -1) fWP_ = po::GetOption<double>("templates.fWP");
  if(fUN_ == -1) fUN_ = po::GetOption<double>("templates.fUN");
};

std::vector<double> IdeogramCombLikelihood::readParameters(const char *whichParameter){
  std::vector<double> pars;
  std::vector<std::string> vsPars;
  std::string sPars = po::GetOption<std::string>(whichParameter);
  boost::split(vsPars, sPars, boost::is_any_of("|"));
  for(auto var : vsPars)
    pars.push_back(std::atof(var.c_str()));
  return pars;
}

void IdeogramCombLikelihood::SetFixedParams(double p0, double p1, double p2, double p3, double p4, double p5, double p6) {
  fp_.push_back(p0);
  fp_.push_back(p1);
  fp_.push_back(p2);
  fp_.push_back(p3);
  fp_.push_back(p4);
  fp_.push_back(p5);
  fp_.push_back(p6);

  useFixedParams_ = true;
}

double IdeogramCombLikelihood::GetFixedParam(int index) {
  return fp_[index];
}

void IdeogramCombLikelihood::SetActive(bool active) {
  active_ = active;
  if (!active) fp_.clear();
}

bool IdeogramCombLikelihood::IsActive() {
  return active_;
}
