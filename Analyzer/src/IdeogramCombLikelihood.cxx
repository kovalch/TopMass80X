#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  
  //* worst case, improvable by f(w_i)
  double fCP = 0.40;
  double fWP = 0.25;
  double fMJ = 0.35;
  //*/

  return p[2] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fMJ * PMJ(x, p));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  double xx = x[0];
  
  double p0 = 1.00528e+01 + 9.34552e-01 * xx;
  
  return TMath::Voigt(p[0] - p0, p[1], 2);
}


namespace cb {
  double A(double alpha, double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(double alpha, double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}


double IdeogramCombLikelihood::PWP(double* x, double* p)
{
  double xx     = x[0];
  
  double N      =  1./0.01;
  double mu     =  8.25340e+00 + 9.02973e-01 * xx;
  double sigma  = -2.77001e+01 + 2.97103e-01 * xx;
  double alpha  = -5.46152e-01 + 5.40408e-03 * xx;
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


double IdeogramCombLikelihood::PMJ(double* x, double* p)
{
  double xx     = x[0];
  
  double N      =  1./0.01;
  double mu     =  4.20677e+01 + 7.22396e-01 * xx;
  double sigma  = -5.18172e+00 + 1.57942e-01 * xx;
  double alpha  = -1.22146e-01 + 3.89904e-03 * xx;
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
