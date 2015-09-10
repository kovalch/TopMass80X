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

    static std::vector<double> muo_parsCP_;
    static std::vector<double> muo_parsWP_;
    static std::vector<double> muo_parsUN_;
    static std::vector<double> muo_parsCPJES_;
    static std::vector<double> muo_parsWPJES_;
    static std::vector<double> muo_parsUNJES_;
    static std::vector<double> ele_parsCP_;
    static std::vector<double> ele_parsWP_;
    static std::vector<double> ele_parsUN_;
    static std::vector<double> ele_parsCPJES_;
    static std::vector<double> ele_parsWPJES_;
    static std::vector<double> ele_parsUNJES_;

    static std::vector<double> muo_massOffset_;
    static std::vector<double> muo_massSlopeMass_;
    static std::vector<double> muo_massSlopeJES_;
    static std::vector<double> muo_massSlopeMassJES_;
    static std::vector<double> muo_jesOffset_;
    static std::vector<double> muo_jesSlopeMass_;
    static std::vector<double> muo_jesSlopeJES_;
    static std::vector<double> muo_jesSlopeMassJES_;
    static std::vector<double> ele_massOffset_;
    static std::vector<double> ele_massSlopeMass_;
    static std::vector<double> ele_massSlopeJES_;
    static std::vector<double> ele_massSlopeMassJES_;
    static std::vector<double> ele_jesOffset_;
    static std::vector<double> ele_jesSlopeMass_;
    static std::vector<double> ele_jesSlopeJES_;
    static std::vector<double> ele_jesSlopeMassJES_;

    static double muo_fCP_, muo_fWP_, muo_fUN_;
    static double ele_fCP_, ele_fWP_, ele_fUN_;
    
    typedef std::tuple<int,double,double> ScanPoint;
    typedef std::map<ScanPoint,double> ScanPointMap; 
    static ScanPointMap PWPnormalizations_;
    static ScanPointMap PUNnormalizations_;


  
  public:
    IdeogramCombLikelihoodLeptonJets();
    ~IdeogramCombLikelihoodLeptonJets() {}
    
    double Evaluate(double *x, double *p);
};

#endif
