#include "TMath.h"

class IdeogramCombLikelihood {
  private:
    double Signal(double *x, double *p);
  
  public:
    IdeogramCombLikelihood() {};
    
    double Evaluate(double *x, double *p);
    double CombBackground(double *x, double *p);
    double CrystalBall(double *x, double *p);
};
