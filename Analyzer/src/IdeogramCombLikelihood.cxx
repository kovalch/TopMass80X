#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool Spring11 = false;
  bool useCalib = false;
  
  //* worst case, improvable by f(w_i)
  double fCP = 0.441627218;
  double fWP = 0.216458497;
  double fUN = 0.341913779;
  //*/
  
  if (onlyCP) {
    fWP = 0;
    fUN = 0;
  }
  
  double MassOffset    = -0.1438;
  double MassSlope     = 0.016;
  
  if (useCalib) {
    x[1] = x[1];
    x[0] = x[0] + MassOffset + MassSlope * (x[0]-172.5);
  }
  
  //return p[0] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[0] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));
}


double IdeogramCombLikelihood::PCP(double *x, double *p) {
  //* hadTopMass
  double mu       = 1.70405e+02 + 9.92056e-01 * (x[0]-172.5) + (7.49913e+01 + 6.49044e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = 9.33475e+00 + 7.75056e-02 * (x[0]-172.5) + (8.32606e+00 + 8.20724e-02 * (x[0]-172.5)) * (x[1]-1.);
  
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
  double mu     =  1.73308e+02 + 7.49383e-01 * (x[0]-172.5) + (8.54287e+01 + 7.46879e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma  =  2.83223e+01 + 2.90576e-01 * (x[0]-172.5) + (2.45982e+01 + 3.59995e-01 * (x[0]-172.5)) * (x[1]-1.);
  double alpha  =  4.73189e-01 + 2.59135e-03 * (x[0]-172.5) + (2.35905e-01 + 4.44735e-03 * (x[0]-172.5)) * (x[1]-1.);
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
  double mu     =  1.67963e+02 + 8.88558e-01 * (x[0]-172.5) + (6.99475e+01 + 5.12963e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma  =  1.94761e+01 + 2.12050e-01 * (x[0]-172.5) + (1.28031e+01 + 2.11590e-01 * (x[0]-172.5)) * (x[1]-1.);
  double alpha  =  7.07848e-01 + 6.83509e-03 * (x[0]-172.5) + (2.39630e-01 + -1.83753e-02 * (x[0]-172.5)) * (x[1]-1.);
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
  double N      =  1;
  
  double mu     =  8.19330e+01 + 1.32156e-02 * (x[0]-172.5) + (4.98160e+01 + -5.74705e-02 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.25837e+00 + 4.22011e-03 * (x[0]-172.5) + (1.39652e+01 + 2.21730e-02 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  6.89876e+00 + -4.83107e-03 * (x[0]-172.5) + (-6.31495e+00 + 5.16539e-02 * (x[0]-172.5)) * (x[1]-1.);
    
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
  
  double mu     =  8.19273e+01 + 1.85945e-02 * (x[0]-172.5) + (5.00479e+01 + 2.54104e-02 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.27086e+00 + 4.13471e-03 * (x[0]-172.5) + (1.49531e+01 + -3.77382e-03 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  7.00045e+00 + -1.26969e-02 * (x[0]-172.5) + (-9.94848e+00 + 9.85874e-03 * (x[0]-172.5)) * (x[1]-1.);
    
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
  
  double mu     =  8.16296e+01 + -4.31147e-03 * (x[0]-172.5) + (3.06442e+01 + -2.97603e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 =  5.64817e+00 + -2.59493e-03 * (x[0]-172.5) + (1.18064e+01 + -1.32256e-01 * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 =  8.38621e+00 + 9.86496e-03 * (x[0]-172.5) + (-9.85783e+00 + 1.69717e-01 * (x[0]-172.5)) * (x[1]-1.);
    
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
