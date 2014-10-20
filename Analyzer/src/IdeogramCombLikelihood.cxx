#include "IdeogramCombLikelihood.h"

#include "ProgramOptionsReader.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

IdeogramCombLikelihood::IdeogramCombLikelihood():
  useFixedParams_(false)
{
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

void IdeogramCombLikelihood::SetFixedParams(double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8) {
  fp_.push_back(p0);
  fp_.push_back(p1);
  fp_.push_back(p2);
  fp_.push_back(p3);
  fp_.push_back(p4);
  fp_.push_back(p5);
  fp_.push_back(p6);
  fp_.push_back(p7);
  fp_.push_back(p8);

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
