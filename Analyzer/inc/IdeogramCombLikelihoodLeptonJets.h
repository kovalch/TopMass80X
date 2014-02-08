#ifndef IDEOGRAMCOMBLIKELIHOODLEPTONJETS_H
#define IDEOGRAMCOMBLIKELIHOODLEPTONJETS_H
#include <iostream>
#include <map>
#include <utility>
#include <vector>

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

    std::vector<double> ele_parsCP_;
    std::vector<double> ele_parsWP_;
    std::vector<double> ele_parsUN_;
    std::vector<double> ele_parsCPJES_;
    std::vector<double> ele_parsWPJES_;
    std::vector<double> ele_parsUNJES_;

    std::vector<double> ele_massOffset_;
    std::vector<double> ele_massSlopeMass_;
    std::vector<double> ele_massSlopeJES_;
    std::vector<double> ele_massSlopeMassJES_;
    std::vector<double> ele_jesOffset_;
    std::vector<double> ele_jesSlopeMass_;
    std::vector<double> ele_jesSlopeJES_;
    std::vector<double> ele_jesSlopeMassJES_;

    double ele_fCP_, ele_fWP_, ele_fUN_;
    
    typedef std::pair<double,double> ScanPoint;
    typedef std::map<ScanPoint,double> ScanPointMap; 
    static ScanPointMap PWPnormalizations_;
    static ScanPointMap PUNnormalizations_;


  
  public:
    IdeogramCombLikelihoodLeptonJets();
    ~IdeogramCombLikelihoodLeptonJets() {}
    
    double Evaluate(double *x, double *p);
};

#endif
