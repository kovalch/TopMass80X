#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "../interface/JECorrectorOnFly.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <stdexcept>

namespace ztop {

double JECorrectorOnFly::getJetEnergyCorrectionValue(const double& jetInitialArea, const double& jetInitialEta, const double& jetInitialPt, const double& rho) {
    
        std::vector<JetCorrectorParameters> vJetPar;
        JetCorrectorParameters *jetParL1, *jetParL2, *jetParL3;
        jetParL1 = new JetCorrectorParameters(filePathL1_);
        jetParL2 = new JetCorrectorParameters(filePathL2_);
        jetParL3 = new JetCorrectorParameters(filePathL3_);
        
        //Order must correspond to one provided by JEC  group: "https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC"
        vJetPar.push_back(*jetParL1);
        vJetPar.push_back(*jetParL2);
        vJetPar.push_back(*jetParL3);
        
        FactorizedJetCorrector *jetCorrector = new FactorizedJetCorrector(vJetPar);
        jetCorrector->setJetA(jetInitialArea);
        jetCorrector->setJetEta(jetInitialEta);
        jetCorrector->setJetPt(jetInitialPt);
        jetCorrector->setRho(rho);
            
        correction_ = jetCorrector->getCorrection();
        delete jetParL1;
        delete jetParL2;
        delete jetParL3;
        delete jetCorrector;
        
        return correction_;
}

void JECorrectorOnFly::setJECorrectionFilesPath(const TString* filePathL1, const TString* filePathL2, const TString* filePathL3) {
    
        filePathL1_ = filePathL1->Data();
        filePathL2_ = filePathL2->Data();
        filePathL3_ = filePathL3->Data();
}

}
