#ifndef IDEOGRAMCOMBLIKELIHOODALLJETS_H
#define IDEOGRAMCOMBLIKELIHOODALLJETS_H
#include <map>
#include <utility>

#include "IdeogramCombLikelihood.h"

#include "TMath.h"

class IdeogramCombLikelihoodAllJets : public IdeogramCombLikelihood {
private:
  double PCP(double *x, double *p);
  double PWP(double *x, double *p);
  double PUN(double *x, double *p);
  
  double PJES(double *x, double *p, double *q);
  //double PCPJES(double *x, double *p, double *q);
  //double PWPJES(double *x, double *p, double *q);
  //double PUNJES(double *x, double *p, double *q);
  
  double PCPJES1(double *x, double *p);
  double PWPJES1(double *x, double *p);
  double PUNJES1(double *x, double *p);
  
  //double PCPJES2(double *x, double *p);
  //double PWPJES2(double *x, double *p);
  //double PUNJES2(double *x, double *p);
  
  double PBKG(double *x, double *p);
  double PBKGJES(double *x, double *p, double *q);
  double PBKGJES1(double *x, double *p);
  //double PBKGJES2(double *x, double *p);
  
  double *parsCP_;
  double *parsWP_;
  double *parsUN_;
  double *parsCPJES_;
  double *parsWPJES_;
  double *parsUNJES_;
  double *parsBKG_;
  double *parsBKGJES_;

  double fSig_, fCP_, fWP_, fUN_;

  double PBKGintegral_;
  
  typedef std::pair<double,double> ScanPoint;
  typedef std::map<ScanPoint,double> ScanPointMap; 
  ScanPointMap PWPintegrals_;
  ScanPointMap PUNintegrals_;

  class MyFunctions {
  public:
    static double langau(Double_t *x, Double_t *p) {
      return p[0]*(p[4]*TMath::Landau(x[0],p[1],p[2],true)+(1-p[4])*TMath::Gaus(x[0],p[1],p[3],true));
    }
    static double lanvog(Double_t *x, Double_t *p) {
      return p[0]*(p[5]*TMath::Landau(x[0],p[1],p[3],true)+(1-p[5])*TMath::Gaus(x[0],p[2],p[4],true));
    }
    static double gamma(Double_t *x, Double_t *p) {
      if(x[0] < p[1]) return 0.;
      else return (TMath::GammaDist(x[0],p[0],p[1],p[2]));
    }
    static double landau(Double_t *x, Double_t *p) {
      return (TMath::Landau(x[0],p[0],p[1]));
    }
    static double gammaland(Double_t *x, Double_t *p) {
      if(x[0] < p[2]) return 0.;
      else return (p[0]/p[6]*TMath::GammaDist(x[0],p[1],p[2],p[3])+(1-p[0])/p[7]*TMath::Landau(x[0],p[4],p[5]));
    }
  };

public:
  IdeogramCombLikelihoodAllJets();
  ~IdeogramCombLikelihoodAllJets();

  double Evaluate(double *x, double *p);
};


#endif
