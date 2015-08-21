#include "IdeogramSampleLikelihood.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

double IdeogramSampleLikelihood::DoEval(const double *x) const {
  double pullWidth(po::GetOption<double>("pullWidth"));
  double pullWidthMuo(po::GetOption<double>("muo_pullWidth"));
  double pullWidthEle(po::GetOption<double>("ele_pullWidth"));
  
  double sampleResult  = 0;
  double sampleSumProb = 0;
  double sampleNEvent  = 0;
  for (const auto& event : eventFunctions_) {
    bool eventIsActive = false;
    int flavour = 0;
    double weight = 1;
    double eventResult  = 0;
    double eventSumProb = 0;
    for (const auto& permutation : event) {
      if (permutation->IsActive()) {
        eventIsActive = true;
        double temp[] = {172.5, 1., 1.,0.5};
        temp[0] = x[0];
        temp[1] = x[1];
        temp[2] = x[2];
        temp[3] = x[3];
        double p[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
        eventResult  += permutation->Evaluate(temp, p);
        eventSumProb += permutation->GetFixedParam(0);
        flavour       = permutation->GetFixedParam(3);
        weight        = permutation->GetFixedParam(8);
      }
    }
    if (eventIsActive) {
      double pullWidthFlavour = pullWidth;
      if     (flavour == 13) pullWidthFlavour = pullWidthMuo;
      else if(flavour == 11) pullWidthFlavour = pullWidthEle;
      //eventSumProb = 1.; // TEST
      sampleResult  += weight * -2.*log(eventResult)*eventSumProb / (pullWidthFlavour*pullWidthFlavour);
      sampleSumProb += weight * eventSumProb;
      sampleNEvent  += weight;
    }
  }
  double constraint = 0;
  if (x[4] > 0.) {
    constraint = pow((x[1]-1.)/x[4], 2);
  }
  return sampleResult * sampleNEvent / sampleSumProb + constraint;
}
