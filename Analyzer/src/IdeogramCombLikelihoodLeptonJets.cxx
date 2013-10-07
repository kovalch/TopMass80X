#include "IdeogramCombLikelihoodLeptonJets.h"

IdeogramCombLikelihoodLeptonJets::IdeogramCombLikelihoodLeptonJets() {};

double IdeogramCombLikelihoodLeptonJets::Evaluate(double *x, double *p) {
  bool onlyCP   = false;
  bool useCalib = false;
  
  //* Summer12
  double fCP = 0.43954; // 1
  double fWP = 0.01067 + 0.188647 + 0.0111505; // 2 3 4
  double fUN = 0.00107329 + 0.0276057 + 0.321233; // -2 -1 6
  //*/
  
  //* Fall11 electrons
  if (p[3] == 11) {
    fCP = 0.447001;
    fWP = 0.00821244 + 0.181794 + 0.00886663;
    fUN = 0.000816016 + 0.0249795 + 0.328331;
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
  double q[12] = {171.674, 0.985702, 79.4061, 0.926074, 9.89094, 0.0800034, 8.80668, -0.00487388, 0, 0, 0, 0};
  double e[12] = {171.646, 0.997112, 80.0779, 0.755107, 9.75059, 0.0779803, 8.97104, -0.0323009, 0, 0, 0, 0};
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
  double q[12] = {173.714, 0.917419, 99.2623, 0.271034, 30.0504, 0.398952, 36.4555, 0.652493, 0.40277, 0.00453923, 0.279563, 0.0017833};
  double e[12] = {174.329, 1.05136, 98.5432, 0.484559, 30.2234, 0.536747, 35.2327, -0.462667, 0.396633, 0.00544819, 0.226293, -0.0123674};
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
  double q[12] = {169.865, 0.87623, 81.0056, 0.835761, 18.697, 0.172046, 7.35387, 0.0221603, 0.831049, 0.00782395, 0.0494347, -0.01188};
  double e[12] = {169.802, 0.898192, 82.3509, 1.05716, 18.5581, 0.177879, 8.47535, -0.2692, 0.818154, 0.00683775, 0.199256, -0.00875282};
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
  double q[12] = {83.1469, 0.00568551, 45.4879, -0.335742, 6.15553, 0.0015377, 12.5085, -0.330609, 7.31215, 0.000182096, -2.84696, 0.209041};
  double e[12] = {83.198, 0.0355267, 41.6371, -0.115806, 6.21283, 0.0138433, 9.42051, 0.0569133, 7.3055, -0.0135193, -0.0826961, 0.288039};
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
  double q[12] = {82.8777, -0.0075446, 45.1709, -0.385762, 6.05136, -0.0018059, 13.4178, -0.321851, 7.51867, 0.00904713, -5.16936, 0.215022};
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
  double q[12] = {82.9979, 0.0228183, 27.1973, 0.382442, 6.90406, 0.00258209, 10.4036, 0.229304, 8.95024, -0.0194394, -7.41695, -0.19067};
  double e[12] = {82.6916, 0.0141625, 27.5887, 0.298442, 6.77016, -0.00468755, 10.7515, 0.135094, 9.14786, -0.0158128, -8.99872, -0.0978321};
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
