#include "IdeogramEventLikelihood.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

double IdeogramEventLikelihood::DoEval(const double *x, const double *p) const {
  double pullWidth(po::GetOption<double>("pullWidth"));
  double pullWidthMuo(po::GetOption<double>("muo_pullWidth"));
  double pullWidthEle(po::GetOption<double>("ele_pullWidth"));
  
  double eventResult  = 0;
  double eventSumProb = 0;
  for (const auto& permutation : eventFunction_) {
    if (permutation->IsActive()) {
      double temp[] = {172.5, 1., 1.,0.5};
      temp[0] = x[0];
      temp[1] = x[1];
      temp[2] = x[2];
      temp[3] = x[3];
      double p[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
      eventResult  += permutation->Evaluate(temp, p);
      eventSumProb += permutation->GetFixedParam(0);
      double flavour       = permutation->GetFixedParam(3);
      double weight        = permutation->GetFixedParam(8);
    }
  }
  return eventResult;
}
