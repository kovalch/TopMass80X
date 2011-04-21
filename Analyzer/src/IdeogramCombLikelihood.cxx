#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 0.3;

  return p[1] * fSig * Signal(x, p) + p[1] * (1.-fSig) * CombBackground(x, p);
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  return TMath::Voigt(p[0] - xx + 0.7, 12, 2);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {

  double xx = x[0];
  
  /*
  double p1 =  38.7578 + 0.756609 * xx;
  double p2 = -46.8933 + 0.434875 * xx;
  double p3 = -69.1613 + 1.86281  * xx;
  double p4 = -102.014 + 0.916367 * xx;
  
  return 2./3. * (TMath::Gaus(p[0], p1, p2, 1) + 1./2. * TMath::Gaus(p[0], p3, p4, 1));
  //*/
  
  ///*
  double p0 =  1.20061 - 0.00383139 * xx;
  double p1 =  190.707 - 0.597429 * xx;
  double p2 = -148.767 + 1.43157 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
  //*/
}
