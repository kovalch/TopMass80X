#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "../interface/JECorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <stdexcept>

namespace ztop {

double JECorrector::getJetEnergyCorrectionValue(const double& jetInitialArea, const double& jetInitialEta, const double& jetInitialPt, const double& rho) {

        jetCorrector_->setJetA(jetInitialArea);
        jetCorrector_->setJetEta(jetInitialEta);
        jetCorrector_->setJetPt(jetInitialPt);
        jetCorrector_->setRho(rho);
            
        correction_ = jetCorrector_->getCorrection();
        
        return correction_;
}

void JECorrector::setJECorrectionFilesPath(const TString* filePathL1, const TString* filePathL2, const TString* filePathL3, const TString* filePathL2L3, const bool& isMC) {
            
        std::vector<JetCorrectorParameters> vJetPar;
        JetCorrectorParameters *jetParL1, *jetParL2, *jetParL3, *jetParL2L3;
        jetParL1 = new JetCorrectorParameters(filePathL1->Data());
        jetParL2 = new JetCorrectorParameters(filePathL2->Data());
        jetParL3 = new JetCorrectorParameters(filePathL3->Data());
        if(!isMC) jetParL2L3 = new JetCorrectorParameters(filePathL2L3->Data());
        
        //Order must correspond to one provided by JEC  group: "https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC"
        vJetPar.push_back(*jetParL1);
        vJetPar.push_back(*jetParL2);
        vJetPar.push_back(*jetParL3);
        if(!isMC) vJetPar.push_back(*jetParL2L3);
        
        jetCorrector_ = new FactorizedJetCorrector(vJetPar);
}

}
