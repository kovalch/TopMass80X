#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 0.3;

  return p[1] * fSig * Signal(x, p) + (1-p[1]) * (1.-fSig) * CombBackground(x, p);
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  return TMath::Gaus(p[0], xx+0.8, 12, 1);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {

  double xx = x[0];
  
  double p0 =  1.27437 - 0.00409434 * xx;
  double p1 =  206.205 - 0.67529 * xx;
  double p2 = -187.656 + 1.65207 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
}
