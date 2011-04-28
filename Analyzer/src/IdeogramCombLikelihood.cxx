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
  
  ///* weight = bProb*fitProb
  double p0 =  5.499088e+00 - 2.810320e-02 * xx;
  if (p0 < 0) p0 = 0;
  double p1 =  8.511500e+01 + 4.477732e-01 * xx;
  double p2 =  20.;
  double p4 =  2.545675e+02 - 2.649908e-01 * xx;
  if (p4 < 0) p4 = 0;
  double p5 =  50.;
  //*/
  
  /* weight = sqrt(bProb*fitProb)
  double p0 =  5.28523e+00 - 2.56035e-02 * xx;
  double p1 =  8.88557e+01 + 4.11095e-01 * xx;
  double p2 =  30.;
  double p4 =  1.81622e+02 + 1.55973e-01 * xx;
  double p5 =  50.;
  //*/
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
  
  /*
  double p0 = -2.23054 + 0.0245958 * xx;
  double p1 =  40.3452 + 0.697798 * xx;
  double p2 = -25.9087 + 0.272108 * xx;
  double p4 = -22.7047 + 1.36774  * xx;
  double p5 =  15.6014 + 0.167157 * xx;
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
  //*/
  
  /*
  double p0 =  2.136350e+00 - 9.930847e-03 * xx;
  double p1 =  3.852494e+02 - 1.807005e+00 * xx;
  double p2 = -3.029468e+02 + 2.404515e+00 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
  //*/
}
