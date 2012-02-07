#include "IdeogramCombLikelihood.h"

double IdeogramCombLikelihood::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool Spring11 = false;
  bool useCalib = false;
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
  
  double JESOffset     = -0.002048;
  double MassOffset    = -0.2359;
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
  double q[8] = {170.84, 0.989457, 81.547, 0.704955, 9.82681, 0.0838884, 9.9775, 0.0157692};
  //* hadTopMass
  double mu       = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma    = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  
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
  double q[12] = {173.608, 0.894004, 97.1005, -0.0998495, 28.4438, 0.412437, 32.6795, -0.303578, 0.469723, 0.0054029, 0.301859, -0.0225543};
  
  double N      =  1./0.01;
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
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
  double q[12] = {168.429, 0.886587, 76.829, 0.410293, 20.0262, 0.283559, 14.3808, 0.0718494, 0.712514, 0.00938252, 0.257174, -0.0138846};
  
  double N      =  1./0.01;
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma  = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double alpha  = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
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
  double q[12] = {82.5221, 0.0399054, 54.3144, 0.0392853, 5.69942, 0.0102082, 15.3069, 0.132278, 7.33692, -0.0179312, -6.52142, 0.0764738};
  
  //* W Mass
  //std::cout << p[4] << " ";
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
   
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
  double q[12] = {82.6462, 0.0144639, 49.5497, 0.610081, 5.82214, -0.00653639, 13.7762, 0.377077, 7.29829, -0.0180972, -6.60668, -0.194419};
  
  //* W Mass
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
      
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
  double q[12] = {81.7715, -0.0109833, 30.4207, -0.43322, 5.98524, -0.00907275, 10.8053, -0.234624, 9.08632, 0.02322, -9.62397, 0.311713};

  //* W Mass
  double N      =  1;
  
  double mu     = q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.);
  double sigma1 = q[4] + q[5] * (x[0]-172.5) + (q[6] + q[7] * (x[0]-172.5)) * (x[1]-1.);
  double sigma2 = q[8] + q[9] * (x[0]-172.5) + (q[10] + q[11] * (x[0]-172.5)) * (x[1]-1.);
    
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
