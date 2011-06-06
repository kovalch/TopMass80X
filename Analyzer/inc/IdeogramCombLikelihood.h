#include "TMath.h"

class IdeogramCombLikelihood {
  private:
    double PCP(double *x, double *p);
    double PWP(double *x, double *p);
    double PMJ(double *x, double *p);
  
  public:
    IdeogramCombLikelihood() {};
    
    double Evaluate(double *x, double *p);
};
