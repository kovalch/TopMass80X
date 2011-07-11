#include <iostream>

#include "TMath.h"

class IdeogramCombLikelihood {
  private:
    double PCP(double *x, double *p);
    double PWP(double *x, double *p);
    double PUN(double *x, double *p);

    double PCPJES(double *x, double *p);
    double PWPJES(double *x, double *p);
    double PUNJES(double *x, double *p);
  
  public:
    IdeogramCombLikelihood() {};
    
    double Evaluate(double *x, double *p);
};
