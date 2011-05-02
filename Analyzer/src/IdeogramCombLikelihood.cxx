#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 1./4.; //TODO

  return p[1] * fSig * Signal(x, p) + (1.-p[1]) * (1.-fSig) * CombBackground(x, p);
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  double p0 =  9.59213e+00 + 9.41938e-01 * xx;
  double p1 =  2.06460e+01 - 6.35914e-02 * xx;
  double p2 = -4.26288e+01 + 2.72355e-01 * xx;
  
  return TMath::Voigt(p[0] - p0, 12*1.1, 2); //TODO
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {
  //TODO: compare different pdfs eventwise over wide range

  double xx = x[0];
  
  double p0 =  1;
  double p1 =  3.45407e+01 + 7.03241e-01 * xx;
  double p2 = -2.33037e+01 + 2.41050e-01 * xx;
  
  return TMath::Landau(p[0], p1, p2, 1);
  // TMath::Landau([0], 3.45407e+01 + 7.03241e-01 * x, -2.33037e+01 + 2.41050e-01 * x)
  
  /*
  double p0 = 0.789244 - 0.00208923 * xx;
  double p1 = 75.;
  double p2 = 9.88908 + 0.589359 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
  //*/
  
  /*
  double p0 =  5.155057e+00 - 2.613158e-02 * xx;
  if (p0 < 0) p0 = 0;
  double p1 =  9.400160e+01 + 3.961189e-01 * xx;
  double p2 =  20.;
  double p4 =  2.403948e+02 - 1.835821e-01 * xx;
  if (p4 < 0) p4 = 0;
  double p5 =  50.;
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
  //*/
}
