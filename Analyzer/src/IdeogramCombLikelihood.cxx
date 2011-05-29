#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  double fSig = 0.8;

  return p[2] * (fSig * Signal(x, p) + (1.-fSig) * CrystalBall(x, p));
}

double IdeogramCombLikelihood::Signal(double *x, double *p) {
  double xx = x[0];
  
  double p0 = 1.00528e+01 + 9.34552e-01 * xx;
  
  return TMath::Voigt(p[0] - p0, p[1], 2);
}

double IdeogramCombLikelihood::CombBackground(double *x, double *p) {
  double xx = x[0];

  /* Landau shape  
  double p0 =  1;
  double p1 =  1.90866e+01 + 8.61702e-01 * xx;
  double p2 = -1.87389e+01 + 2.42629e-01 * xx;
  
  return TMath::Landau(p[0], p1, p2, 1);
  //*/
  // TMath::Landau([0], 3.45407e+01 + 7.03241e-01 * x, -2.33037e+01 + 2.41050e-01 * x)
  
  /* LogNormal shape
  double p0 = 0.789244 - 0.00208923 * xx;
  double p1 = 75.;
  double p2 = 9.88908 + 0.589359 * xx;
  
  return TMath::LogNormal(p[0], p0, p1, p2);
  //*/
  
  /* double gaus shape
  double p0 =  4.10310e+00 - 2.01348e-02 * xx;
  if (p0 < 0) p0 = 0;
  double p1 =  5.85993e+01 + 6.02224e-01 * xx;
  double p2 =  20.;
  double p4 =  1.65129e+02 + 2.56381e-01 * xx;
  if (p4 < 0) p4 = 0;
  double p5 =  50.;
  
  return 1./(p0 + 1.) * (p0 * TMath::Gaus(p[0], p1, p2, 1) + TMath::Gaus(p[0], p4, p5, 1));
  //*/
}

namespace cb {
  double A(double alpha, double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(double alpha, double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}

// 7 parameters: [0] -> [6]
double IdeogramCombLikelihood::CrystalBall(double* x, double* p)
{
  double xx     = x[0];
  
  double N      =  1./0.01;
  double mu     =  3.56187e+01 + 7.53376e-01 * xx;
  double sigma  = -2.23761e+01 + 2.60009e-01 * xx;
  double alpha  = -1.40038e-01 + 3.53449e-03 * xx;
  double power  =  15;
  double t = (p[0] - mu) / sigma;
  
  //*
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (-1+TMath::Erf(-sigma*alpha)/(sqrt(2)*sigma));
  double N2 = cb::A(alpha,power) * sigma/(-1.+power) * ( 0. + //lim(x->inf)
                pow(((cb::B(alpha,power) * sigma - sigma*alpha)/sigma), (1. - power))
              );
  N = N1 + N2;
  //*/

  if(t < alpha)
    return 1./N * TMath::Exp(-t*t/2);
  else
    return 1./N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}
