#include "MassAnalyzer.h"

MassAnalyzer::MassAnalyzer(const std::string& identifier, TTree* tree) :
  fIdentifier_(identifier), fTree_(tree)
{
}

std::pair<double, double>
MassAnalyzer::GetValue(const std::string& whichValue) const {
  std::map<std::string, std::pair<double, double> >::const_iterator value_iterator = values_.find(whichValue);
  if(value_iterator != values_.end()){
    return value_iterator->second;
  }
  else{
    std::cerr << "WARNING: The searched value *" << whichValue << "* does not exist!" << std::endl;
    assert(0);
  }
  return std::make_pair(-1, -1.);
}

std::map<std::string, std::pair<double, double> >
MassAnalyzer::GetValues() const{
  return values_;
}

void
MassAnalyzer::SetValue(const std::string& whichValue, double val, double valError){
  values_[whichValue] = std::make_pair(val, valError);
}
