#ifndef IDEOGRAMCOMBLIKELIHOODLEPTONJETS_H
#define IDEOGRAMCOMBLIKELIHOODLEPTONJETS_H
#include <iostream>

#include "IdeogramCombLikelihood.h"

#include "TMath.h"

class IdeogramCombLikelihoodLeptonJets : public IdeogramCombLikelihood {
  private:
    double PCP(double *x, double *p);
    double PWP(double *x, double *p);
    double PUN(double *x, double *p);

    double PCPJES(double *x, double *p);
    double PWPJES(double *x, double *p);
    double PUNJES(double *x, double *p);
  
  public:
    IdeogramCombLikelihoodLeptonJets();
    ~IdeogramCombLikelihoodLeptonJets() {}
    
    double Evaluate(double *x, double *p);
};

#endif
