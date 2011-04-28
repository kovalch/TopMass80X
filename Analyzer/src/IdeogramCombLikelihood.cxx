#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 1./12.;

  return p[1] * (fSig * Signal(x, p) + (1.-fSig) * CombBackground(x, p));
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  return TMath::Voigt(p[0] - xx + 0.7, 12*1.1, 2);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {

  double xx = x[0];
  
  double p0 =  5.155057e+00 - 2.613158e-02 * xx;
  if (p0 < 0) p0 = 0;
  double p1 =  9.400160e+01 + 3.961189e-01 * xx;
  double p2 =  20.;
  double p4 =  2.403948e+02 - 1.835821e-01 * xx;
  if (p4 < 0) p4 = 0;
  double p5 =  50.;
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
}
