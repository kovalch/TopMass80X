#ifndef IDEOGRAMCOMBLIKELIHOODALLJETS_H
#define IDEOGRAMCOMBLIKELIHOODALLJETS_H
#include <map>
#include <utility>
#include <vector>

#include "IdeogramCombLikelihood.h"

#include "TMath.h"

class IdeogramCombLikelihoodAllJets : public IdeogramCombLikelihood {
private:
  double PCP(double *x, double *p);
  double PWP(double *x, double *p);
  double PUN(double *x, double *p);
  
  double PJES(double *x, double *p, const std::vector<double> &q);
  
  double PCPJES1(double *x, double *p);
  double PWPJES1(double *x, double *p);
  double PUNJES1(double *x, double *p);
  
  double PBKG(double *x, double *p);
  double PBKGJES(double *x, double *p, const std::vector<double> &q);
  double PBKGJES1(double *x, double *p);
  
  std::vector<double> readParameters(const char *whichParameter);

  static std::vector<double> parsCP_;
  static std::vector<double> parsWP_;
  static std::vector<double> parsUN_;
  static std::vector<double> parsCPJES_;
  static std::vector<double> parsWPJES_;
  static std::vector<double> parsUNJES_;
  static std::vector<double> parsBKG_;
  static std::vector<double> parsBKGJES_;

  static std::vector<double> massOffset_;
  static std::vector<double> massSlopeMass_;
  static std::vector<double> massSlopeJES_;
  static std::vector<double> massSlopeMassJES_;
  static std::vector<double> jesOffset_;
  static std::vector<double> jesSlopeMass_;
  static std::vector<double> jesSlopeJES_;
  static std::vector<double> jesSlopeMassJES_;

  static double fSig_, fCP_, fWP_, fUN_;

  static double PBKGintegral_;
  
  typedef std::pair<double,double> ScanPoint;
  typedef std::map<ScanPoint,double> ScanPointMap; 
  static ScanPointMap PWPnormalizations_;
  static ScanPointMap PUNnormalizations_;

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

public:
  IdeogramCombLikelihoodAllJets();
  ~IdeogramCombLikelihoodAllJets();

  double Evaluate(double *x, double *p);
};


#endif
