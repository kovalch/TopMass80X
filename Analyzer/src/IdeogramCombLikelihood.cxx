#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 0.;

  return p[1] * (fSig * Signal(x, p) + (1.-fSig) * CombBackground(x, p));
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  return TMath::Voigt(p[0] - xx + 0.7, 12, 2);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {

  double xx = x[0];
  
  ///*
  double p0 = -2.23054 + 0.0245958 * xx;
  double p1 =  40.3452 + 0.697798 * xx;
  double p2 = -25.9087 + 0.272108 * xx;
  double p4 = -22.7047 + 1.36774  * xx;
  double p5 =  15.6014 + 0.167157 * xx;
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
  //*/
  
  /*
  double p0 =  1.20061 - 0.00383139 * xx;
  double p1 =  190.707 - 0.597429 * xx;
  double p2 = -148.767 + 1.43157 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
  //*/
}
