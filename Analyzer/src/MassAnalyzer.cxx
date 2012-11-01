#include "MassAnalyzer.h"

MassAnalyzer::MassAnalyzer(const TString& identifier, TTree* tree) : 
  fIdentifier_(identifier), fTree_(tree)
{
}

std::pair<double, double>
MassAnalyzer::GetValue(TString whichValue) const {
  std::map<TString, std::pair<double, double> >::const_iterator value_iterator = values_.find(whichValue);
  if(value_iterator != values_.end()){
    return value_iterator->second;
  }
  else{
    std::cerr << "WARNING: The searched value *" << whichValue << "* does not exist!" << std::endl;
    assert(0);
  }
  return std::make_pair(-1, -1.);
}

std::map<TString, std::pair<double, double> >
MassAnalyzer::GetValues() const{
  return values_;
}

void
MassAnalyzer::SetValue(TString whichValue, double val, double valError){
  values_[whichValue] = std::make_pair(val, valError);
}
