#ifndef correctJE_h
#define correctJE_h

#include "TString.h"

/*
 WHATEVER you add , please don't use exit() in case an error occurs.
 replace it with either:
 - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 - return something (-1 or another int for debugging)
*/

class FactorizedJetCorrector;

namespace ztop {

class JECorrector {

public:
    JECorrector():correction_(1.0){}
    ~JECorrector(){}

    double getJetEnergyCorrectionValue(const double& jetInitialArea, const double& jetInitialEta, const double& jetInitialPt, const double& rho);
    void setJECorrectionFilesPath(const TString* filePathL1, const TString* filePathL2, const TString* filePathL3, const TString* filePathL2L3, const bool& isMC);

private:
    double correction_;
    FactorizedJetCorrector *jetCorrector_;
    
};

}

#endif
