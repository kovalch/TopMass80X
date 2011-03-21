#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 0.7;
  double xx = x[0];

  return p[2] * fSig * Signal(x, p) + (1-fSig) * CombBackground(x, p);
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  return TMath::Gaus(p[0],xx,p[1],1);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {
  return 1/(3.5+1) * (3.5 * TMath::Gaus(p[0], (6.329e+01 + p[0]*0.6317), 25, 1) + TMath::Gaus(p[0], 225, 46, 1));
}
