#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  
  //* worst case, improvable by f(w_i)
  double fCP = 0.404;
  double fWP = 0.246;
  double fUN = 0.350;
  //*/
  
  fWP = 0;
  fUN = 0;
  
  //std::cout << "x[0]: " << x[0] << " x[1]: " << x[1] << std::endl;
  //std::cout << p[2] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p)) << std::endl;
  
  //return p[2] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[2] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  double mu       = 1.71262e+02 + 9.44736e-01 * (x[0]-172.5) + (x[1]-1.) * (7.84333e+01 + 0*9.31738e-01 * (x[0]-172.5));
  //double mu = 1.60174e+01 + 8.93327e-01 * x[0] + (x[1]-1.) * (7.84333e+01 + 9.31738e-01 * (x[0]-172.5));
  double sigma    = p[1]                                    ;// + (x[1]-1.) * (8.81344e+00);// - 1.56275e-01 * (x[0]-172.5));
  
  //return 1; // JES test
  return TMath::Voigt(p[0] - mu, sigma, 2);
  //return TMath::Gaus(p1, p[3]/80.4, 0.1);
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
  double mu     =  2.94263e+01 + 8.05300e-01 * x[0] + (x[1]-1.) * 1.11775e+02;
  double sigma  = -2.49181e+01 + 2.81900e-01 * x[0] + (x[1]-1.) * 3.74308e+01;
  double alpha  =  5.79839e-04 + 2.40480e-03 * x[0] + (x[1]-1.) * 4.14620e-01;
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
  double mu     =  5.27110e+01 + 6.71607e-01 * x[0] + (x[1]-1.) * 8.19469e+01;
  double sigma  =  1.41279e+00 + 1.15346e-01 * x[0] + (x[1]-1.) * 1.99421e+01;
  double alpha  = -1.95741e-01 + 5.00641e-03 * x[0] + (x[1]-1.) * 2.03677e-01;
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
  /* W Mass
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
  
  //* pt balance
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
  
  double mu     =  8.18935e+01 + 3.88707e+01 * (x[1]-1.);
  double sigma1 =  5.45248e+00 + 8.17274e+00 * (x[1]-1.);
  double sigma2 =  7.29203e+00 - 5.61393e+00 * (x[1]-1.);
    
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
  
  double mu     =  8.15940e+01 + 1.81073e+01 * (x[1]-1.);
  double sigma1 =  5.52606e+00 + 2.35925e+00 * (x[1]-1.);
  double sigma2 =  8.36194e+00 - 8.36194e+00 * (x[1]-1.);
    
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
