#include "IdeogramSampleLikelihood.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

double IdeogramSampleLikelihood::DoEval(const double *x) const {
  double pullWidth(po::GetOption<double>("pullWidth"));
  
  double sampleResult  = 0;
  double sampleSumProb = 0;
  double sampleNEvent  = 0;
  for (const auto& event : eventFunctions_) {
    double eventResult  = 0;
    double eventSumProb = 0;
    for (const auto& permutation : event) {
      double temp[] = {172.5, 1.};
      temp[0] = x[0];
      temp[1] = x[1];
      double p[] = {0., 0., 0., 0., 0., 0., 0.};
      eventResult  += permutation->Evaluate(temp, p);
      eventSumProb += permutation->GetFixedParam(0);
    }
    sampleResult  += -2.*log(eventResult)*eventSumProb;
    sampleSumProb += eventSumProb;
    sampleNEvent  += 1.;
  }
  return sampleResult * sampleNEvent / sampleSumProb / (pullWidth*pullWidth);
}

