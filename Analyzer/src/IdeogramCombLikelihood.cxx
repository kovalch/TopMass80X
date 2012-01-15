#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool Spring11 = false;
  bool useCalib = true;
  bool useFix   = false;
  
  //* worst case, improvable by f(w_i)
  double fCP = 0.444401948;
  double fWP = 0.21610251;
  double fUN = 0.339498595;
  //*/
  
  if (onlyCP) {
    fWP = 0;
    fUN = 0;
  }
  
  double JESOffset     = -0.002026;
  double MassOffset    = -0.3623;
  double MassSlope     = 0.016;
  
  double JESFix  = 0.00333;
  double MassFix = -0.249;
  
  if (useCalib) {
    x[1] = x[1] + JESOffset;
    x[0] = x[0] + MassOffset; // + MassSlope * (x[0]-172.5);
  }
  
  if (useFix) {
    x[1] = x[1] + JESFix;
    x[0] = x[0] + MassFix; // + MassSlope * (x[0]-172.5);
  }
  
  //return p[0] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[0] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));// * TMath::Gaus(x[1], 1, 0.01, kTRUE);
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  //* hadTopMass
  double mu       = 1.70821e+02 + 9.91801e-01 * (x[0]-172.5) + (8.15982e+01 + 6.84800e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = 9.75238e+00 + 8.29301e-02 * (x[0]-172.5) + (9.65337e+00 + 3.52131e-02 * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[1] - mu, sigma, 2);
  //*/
  
  /* sin(theta*)
  double mu       = 1.13011e+00 - 1.15254e-02 * (x[0]-172.5) - 5.50104e-01 * (x[1]-1.);
  double sigma    = 0.3;
  
  return TMath::Landau(p[1], mu, sigma, true);
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


double IdeogramCombLikelihood::PWP(double* x, double* p)
{
  double N      =  1./0.01;
  double mu     =  1.77356e+02 + 1.03758e+00 * (x[0]-172.5) + (9.51202e+01 - 6.67208e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma  =  2.90139e+01 + 4.20973e-01 * (x[0]-172.5) + (2.83814e+01 - 8.44371e-01 * (x[0]-172.5)) * (x[1]-1.);
  double alpha  =  4.77378e-01 + 5.26067e-03 * (x[0]-172.5) + (2.03662e-01 - 3.60043e-02 * (x[0]-172.5)) * (x[1]-1.);
  double power  =  15;
  double t = (p[1] - mu) / sigma;
  
  //*
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
  double N2 = cb::A(alpha,power) / (-1.+power) * (
                (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
               -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
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
  double mu     =  1.68434e+02 + 8.99773e-01 * (x[0]-172.5) + (7.70014e+01 + 3.51422e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma  =  1.97809e+01 + 2.77859e-01 * (x[0]-172.5) + (1.42499e+01 + 1.39767e-01 * (x[0]-172.5)) * (x[1]-1.);
  double alpha  =  7.04974e-01 + 9.32045e-03 * (x[0]-172.5) + (2.85543e-01 - 1.64546e-02 * (x[0]-172.5)) * (x[1]-1.);
  double power  =  5;
  double t = (p[1] - mu) / sigma;
  
  //*
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
  double N2 = cb::A(alpha,power) / (-1.+power) * (
                (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
               -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
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
  //std::cout << p[4] << " ";
  double N      =  1;
  
  double mu     =  8.23959e+01 + 4.59395e-02 * (x[0]-172.5) + (5.12109e+01 + 5.77990e-02 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.41708e+00 + 1.57805e-02 * (x[0]-172.5) + (1.44012e+01 + 1.43676e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  7.03401e+00 - 1.99430e-02 * (x[0]-172.5) + (-6.85933e+00 + 7.15651e-02 * (x[0]-172.5)) * (x[1]-1.);
    
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
  double mu       = 1.36597e-02 - 1.00417e-01 * (x[1]-1.);
  double sigma    = 4.41279e-02 + 3.22201e-02 * (x[1]-1.);
  double width    = 6.67793e-02 - 4.12224e-02 * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
	
	/* Constraint
	return TMath::Gaus(x[1], 1, 0.1, kTRUE);
	//*/
	
	/* b scale estimator
  double mu       = 1.23115e+00 - 4.21823e-03 * (x[0]-172.5) + (-1.48697e+00 - 4.92782e-03 * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = 0.001;
  double width    = 1.29952e+00 + 1.10928e-02 * (x[0]-172.5) + (-1.14917e+00 - 3.74889e-02 * (x[0]-172.5)) * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}

double IdeogramCombLikelihood::PWPJES(double* x, double* p)
{
  //* W Mass
  double N      =  1;
  
  double mu     =  8.29371e+01 + 4.36646e-03 * (x[0]-172.5) + (4.17956e+01 + 4.02124e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.76377e+00 - 1.28977e-02 * (x[0]-172.5) + (9.09588e+00 + 2.28455e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  6.69985e+00 - 1.11469e-02 * (x[0]-172.5) + (-3.44082e+00 - 2.47573e-03 * (x[0]-172.5)) * (x[1]-1.);
    
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
  double mu       = 1.94840e-02 - 5.30632e-02 * (x[1]-1.);
  double sigma    = 5.86415e-02 + 2.58850e-02 * (x[1]-1.);
  double width    = 1.22249e-01 - 1.36272e-01 * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
  
  /* b scale estimator
  double mu       = 7.18283e-02 + 1.68825e-01 * (x[1]-1.);
  double sigma    = 9.92390e-01 - 8.84443e-01 * (x[1]-1.);
  double width    = 2.15202e+00 - 2.33270e+00 * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}

double IdeogramCombLikelihood::PUNJES(double* x, double* p)
{
  //* W Mass
  double N      =  1;
  
  double mu     =  8.18050e+01 - 1.70440e-02 * (x[0]-172.5) + (2.84184e+01 - 3.79385e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.70154e+00 - 1.29138e-02 * (x[0]-172.5) + (1.01906e+01 - 2.63053e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  8.45273e+00 + 2.61732e-02 * (x[0]-172.5) + (-8.36431e+00 + 2.19731e-01 * (x[0]-172.5)) * (x[1]-1.);
    
  double t1     =  (p[2] - mu) / sigma1;
  double t2     =  (p[2] - mu) / sigma2;
  
  double N1     =  sqrt(2*sigma1*sigma1);
  double N2     =  sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2.;
  
  if(t1 < 0)
    return 1./N * TMath::Exp(-t1*t1/2);
  else
    return 1./N * TMath::Exp(-t2*t2/2);
  //*/
  
  /* pt balance
  double mu       = 9.99882e-03 + 2.33525e-02 * (x[1]-1.);
  double sigma    = 5.51757e-02 - 2.62315e-02 * (x[1]-1.);
  double width    = 1.35272e-01 - 3.44117e-03 * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
  
  /* b scale estimator
  double mu       = 9.45273e-01 - 1.07859e+00 * (x[1]-1.) - 0.2;
  double sigma    = 0.001;
  double width    = 2.14746e+00 - 2.29404e+00 * (x[1]-1.);
  
  return TMath::Voigt(p[3] - mu, sigma, width);
  //*/
}
