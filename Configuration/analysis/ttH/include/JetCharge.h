#ifndef JetCharge_h
#define JetCharge_h

#include <vector>
#include <map>

class TString;

#include "AnalyzerBase.h"
#include "../../common/include/classesFwd.h"
#include "../../common/include/storeTemplate.h"

class JetCategories;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinematicReconstructionSolutions;
namespace tth{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}


/// jet track reweighting class, test here whatever you want to test
class JetCharge {
    
public:
    
    /// Constructor
    JetCharge();
    
    /// Destructor
    ~JetCharge(){}
    
    /// p_{T} weighted jet charge calculation for a given squeezing parameter x (optimal value x = 0.8)
    double ptWeightedJetChargeX (const int jetId, const LV& recoJet, const double x, const std::vector<int> pfCandidateJetIndex, const VLV& pfCandidates, const std::vector<int> pfCandidateCharge);
    
};

#endif