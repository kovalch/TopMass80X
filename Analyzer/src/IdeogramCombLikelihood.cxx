#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  
  //* worst case, improvable by f(w_i)
  double fCP = 0.390;
  double fWP = 0.235;
  double fUN = 0.375;
  //*/
  
  //std::cout << "x[0]: " << x[0] << " x[1]: " << x[1] << std::endl;
  //std::cout << p[2] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p)) << std::endl;
  
  //return p[2] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[2] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  double mu       = 1.72119e+02 + 9.36446e-01 * (x[0]-172.5) + 7.91540e+01 * (x[1]-1.);
  double sigma    = p[1];
  
  return TMath::Voigt(p[0] - mu, sigma, 2);
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
  double N      =  1./0.01;
  double mu     =  1.68182e+02 + 9.81981e-01 * (x[0]-172.5) + 1.10020e+02 * (x[1]-1.);
  double sigma  =  2.41622e+01 + 4.69440e-01 * (x[0]-172.5) + 3.54726e+01 * (x[1]-1.);
  double alpha  =  4.34182e-01 + 4.35156e-03 * (x[0]-172.5) + 4.09760e-01 * (x[1]-1.);
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


double IdeogramCombLikelihood::PUN(double* x, double* p)
{
  double xx     = x[0];
  
  double N      =  1./0.01;
  double mu     =  1.68631e+02 + 7.65029e-01 * (x[0]-172.5) + 7.77693e+01 * (x[1]-1.);
  double sigma  =  2.18236e+01 + 2.38588e-01 * (x[0]-172.5) + 2.76385e+01 * (x[1]-1.);
  double alpha  =  7.09159e-01 + 7.24736e-03 * (x[0]-172.5) + 3.67483e-01 * (x[1]-1.);
  double power  =  5;
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

double IdeogramCombLikelihood::PCPJES(double* x, double* p)
{
  //* W Mass
  double N      =  1;
  
  double mu     =  8.25053e+01 + 4.63537e+01 * (x[1]-1.);
  double sigma1 =  5.65204e+00 + 1.15149e+01 * (x[1]-1.);
  double sigma2 =  7.21349e+00 - 2.84521e+00 * (x[1]-1.);
    
  double t1     =  (p[3] - mu) / sigma1;
  double t2     =  (p[3] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
  double N      =  1./0.01;
  double mu     =  3.36513e+00 + 5.14285e+01 * (x[1]-1.);
  double sigma  =  1.74507e+01 + 1.75185e+01 * (x[1]-1.);
  double alpha  =  1.38928e+00;
  double power  =  3;
  double t = (p[3] - mu) / sigma;
  
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (-1+TMath::Erf(-sigma*alpha)/(sqrt(2)*sigma));
  double N2 = cb::A(alpha,power) * sigma/(-1.+power) * ( 0. + //lim(x->inf)
                pow(((cb::B(alpha,power) * sigma - sigma*alpha)/sigma), (1. - power))
              );
  N = N1 + N2;
  
  if(t < alpha)
    return 1./N * TMath::Exp(-t*t/2);
  else
    return 1./N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
  //*/
}

double IdeogramCombLikelihood::PWPJES(double* x, double* p)
{
  double N      =  1;
  
  double mu     =  8.24119e+01 + 3.67810e+01 * (x[1]-1.);
  double sigma1 =  5.55929e+00 + 7.91431e+00 * (x[1]-1.);
  double sigma2 =  7.49184e+00 - 5.54074e+00 * (x[1]-1.);
    
  double t1     =  (p[3] - mu) / sigma1;
  double t2     =  (p[3] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
}

double IdeogramCombLikelihood::PUNJES(double* x, double* p)
{
  double N      =  1;
  
  double mu     =  8.20562e+01 + 2.18725e+01 * (x[1]-1.);
  double sigma1 =  5.72921e+00 + 7.09411e+00 * (x[1]-1.);
  double sigma2 =  8.34403e+00 - 6.50901e+00 * (x[1]-1.);
    
  double t1     =  (p[3] - mu) / sigma1;
  double t2     =  (p[3] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
}
