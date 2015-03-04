#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <map>
#include <vector>

#include <TString.h>
#include <Math/VectorUtil.h>
#include <TProfile.h>
#include <TObjString.h>

#include "JetCharge.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "higgsUtils.h"
#include "HiggsAnalysis.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/storeTemplate.h"
#include "../../common/include/classesFwd.h"


JetCharge::JetCharge() {}


double JetCharge::ptWeightedJetChargeX (const int jetId, const LV& recoJet, const double x, const std::vector<int> pfCandidateJetIndex, const VLV& pfCandidates, const std::vector<int> pfCandidateCharge)
{
    // Access jet momentum information
    double jetPx = recoJet.px();
    double jetPy = recoJet.py();
    double jetPz = recoJet.pz();
    
    // Define relevant variables for c_{rel} calculation
    double sumMomentum = 0.;
    double sumMomentumQ = 0.;
    
    for (size_t iCandidate=0;iCandidate!=pfCandidates.size();++iCandidate)
    {
        // Check that the pfCandidate corresponds to the jet
        if (jetId!=pfCandidateJetIndex.at(iCandidate)) continue;
        
        // Access pfCandidate mometum and charge information
        const double constituentPx = pfCandidates.at(iCandidate).px();
        const double constituentPy = pfCandidates.at(iCandidate).py();
        const double constituentPz = pfCandidates.at(iCandidate).pz();
        const double product = constituentPx*jetPx + constituentPy*jetPy + constituentPz*jetPz;
        
        int charge = pfCandidateCharge.at(iCandidate);
        
        // Sum over all the pfCandidates
        const double productPow = std::pow(product, x);
        sumMomentum += productPow;
        sumMomentumQ += charge*productPow;
    }
    
    // Obtain the jet c_{rel}
    const double ptWeightedJetChargeXValue(sumMomentum>0 ? sumMomentumQ/sumMomentum : 0);
    return ptWeightedJetChargeXValue;
}
    