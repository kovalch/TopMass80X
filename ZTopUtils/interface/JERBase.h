#ifndef JERBASE_H
#define JERBASE_H


#include <vector>
#include "../interface/miscUtils.h"
#include <algorithm>
#include <iostream>

// following recommendation https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution;   https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefSyst

namespace ztop {

// typedef ROOT::Math::PolarLorentzVector<ROOT::Math::PxPyPzE4D<double> > PolarLorentzVector;
/**
 *
 WHATEVER you add as functions, please don't use exit() in case an error occurs.
 replace it with either:
 - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 - return something (-1 or another int for dubugging)

 */

class JERBase {
public:
    JERBase() {
        setSystematics("def");
    }
    ~JERBase() {
    }

    void setResolutionEtaRanges(std::vector<float> ranges) {
        resranges_ = ranges;
    }
    void setResolutionFactors(std::vector<float> factors) {
        resfactors_ = factors;
    }

    void setSystematics(std::string type);

    void correctP4(float & recopt, float& recoabs_eta, float & recophi, float & recom, //full lorentzvector
            const float & genpt) const;

protected:

    std::vector<float> resranges_;
    std::vector<float> resfactors_;

};
}
#endif
