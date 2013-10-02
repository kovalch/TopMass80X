#include "IdeogramCombLikelihoodLeptonJets.h"

IdeogramCombLikelihoodLeptonJets::IdeogramCombLikelihoodLeptonJets() {};

double IdeogramCombLikelihoodLeptonJets::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool Spring11 = false;
  bool useCalib = true;
  
  //* Fall11
  double fCP = 0.438236;
  double fWP = 0.207301;
  double fUN = 0.354463;
  //*/
  
  //* Fall11 electrons
  if (p[3] == 11) {
    fCP = 0.439947;
    fWP = 0.205682;
    fUN = 0.354370;
  }
  //*
  
  /* Fall11 P11
  double fCP = 0.415661;
  double fWP = 0.201059;
  double fUN = 0.38328;
  //*/
  
  /* worst case, improvable by f(w_i)
  double fCP = 0.439252;
  double fWP = 0.211647;
  double fUN = 0.349102;
  //*/
  
  if (onlyCP) {
    fWP = 0;
    fUN = 0;
  }
  
  double MassOffset       = -9.58404e-02;
  double MassSlopeMass    = -1.82737e-02;
  double MassSlopeJES     =  5.25560e+00;
  double MassSlopeMassJES = -7.13269e-02;
                            
  double JESOffset        = -3.06117e-03;
  double JESSlopeMass     =  2.87550e-04;
  double JESSlopeJES      = -3.89880e-02;
  double JESSlopeMassJES  =  1.19222e-03;
  
  if (p[3] == 11) {
    MassOffset       = -3.60599e-01;
    MassSlopeMass    = -2.26131e-03;
    MassSlopeJES     = -2.66705e+00;
    MassSlopeMassJES = -1.84892e-01;
                       
    JESOffset        = -1.58004e-03;
    JESSlopeMass     = -7.27182e-05;
    JESSlopeJES      =  1.24220e-02;
    JESSlopeMassJES  =  8.35960e-04;
  }
  
  double m = x[0];
  double j = x[1];
  
  if (useCalib) {
    x[1] = j + JESOffset + JESSlopeMass*(m-172.5) + JESSlopeJES*(j-1.) + JESSlopeMassJES*(m-172.5)*(j-1.);
    x[0] = m + MassOffset + MassSlopeMass*(m-172.5) + MassSlopeJES*(j-1.) + MassSlopeMassJES*(m-172.5)*(j-1.);
  }
  
  //return p[0] * (fCP * PCP(x, p) + fWP * PWP(x, p) + fUN * PUN(x, p));
  return p[0] * (fCP * PCP(x, p) * PCPJES(x, p) + fWP * PWP(x, p) * PWPJES(x, p) + fUN * PUN(x, p) * PUNJES(x, p));
}


double IdeogramCombLikelihoodLeptonJets::PCP(double *x, double *p) {
  double q[12] = {170.729, 0.980962, 86.9952, 0.810794, 10.193, 0.0739653, 10.589, -0.103049, 0, 0, 0, 0};
  double e[12] = {170.772, 0.950394, 87.4086, 0.776666, 10.1953, 0.0950286, 11.012, 0.195389, 0, 0, 0, 0};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];
  
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


double IdeogramCombLikelihoodLeptonJets::PWP(double* x, double* p)
{
  double q[12] = {172.06, 0.720447, 90.3707, -0.860654, 27.934, 0.311025, 29.8054, -1.07054, 0.426473, 0.00193419, 0.125964, 0.00822953};
  double e[12] = {174.123, 0.620424, 96.8582, 0.885688, 29.0981, 0.140688, 30.3437, 0.734781, 0.435573, -0.0022993, 0.123777, -0.0108803};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];
  
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


double IdeogramCombLikelihoodLeptonJets::PUN(double* x, double* p)
{
  double q[12] = {168.581, 0.887245, 87.4702, 0.645641, 20.4999, 0.249314, 18.5896, 0.155375, 0.876613, 0.00951753, 0.319, 0.00182893};
  double e[12] = {168.562, 0.842331, 86.132, 0.77473, 20.4059, 0.250519, 18.0985, 0.442554, 0.864017, 0.00832897, 0.218885, 0.0219663};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];
  
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

double IdeogramCombLikelihoodLeptonJets::PCPJES(double* x, double* p)
{
  double q[12] = {82.7622, 0.0296313, 54.7387, 0.261219, 5.7984, 0.002663, 15.5537, 0.186603, 7.23293, -0.0154696, -7.29268, -0.0437143};
  double e[12] = {82.8179, -0.0130562, 54.298, 0.26241, 5.79095, -0.0144692, 15.1857, 0.157664, 7.23174, 0.0149186, -6.9936, -0.240175};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];
  
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

double IdeogramCombLikelihoodLeptonJets::PWPJES(double* x, double* p)
{
  double q[12] = {82.5335, 0.0134514, 52.2913, 0.215036, 5.78653, -0.0171043, 15.1453, 0.010765, 7.46176, -0.0201443, -9.53141, 0.100091};
  double e[12] = {82.5788, 0.0143882, 53.0209, 0.00501586, 5.73298, 0.00455929, 15.6096, -0.0153593, 7.43306, -0.0127021, -9.1285, 0.0178314};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];
  
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

double IdeogramCombLikelihoodLeptonJets::PUNJES(double* x, double* p)
{
  double q[12] = {82.2989, 0.0220154, 32.6116, -0.218329, 6.23978, 0.00696407, 11.4367, -0.0979023, 8.74599, -0.00031112, -10.8387, 0.0516436};
  double e[12] = {82.3831, 0.00858293, 31.3081, 1.05433, 6.274, -0.00778027, 10.5673, 0.575389, 8.77308, -0.0075097, -9.95694, -0.571045};
  if (p[3] == 11) for (int i = 0; i<12; ++i) q[i] = e[i];

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
