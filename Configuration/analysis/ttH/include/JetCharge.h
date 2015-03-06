#ifndef JetCharge_h
#define JetCharge_h

#include <vector>
#include <map>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
    
    /// p weighted jet charge calculation for a given squeezing parameter x (optimal value x = 0.8)
    double ptWeightedJetChargeX (const int jetId, const LV& recoJet, const double x, const std::vector<int> pfCandidateJetIndex, const VLV& pfCandidates, const std::vector<int> pfCandidateCharge);
   
    /// MVA charge weight -> This weight should be understood as a charge by the top-system identification mva
    double mvaJetChargeWeight(const int jetIndex, const LV& jet, const RecoObjects& recoObjects);
    
    TMVA::Reader* jetChargeReader_;
    TString case1_;
    
    TTree* mvaJetChargeTestTree_;
    TTree* mvaJetChargeTrainTree_;
    
    struct MvaJetChargeVariables
    {
        Float_t longChargeJet_;
        Float_t relChargeJet_;
        Float_t leadingTrackPtWeightedCharge_;
        Float_t subleadingTrackPtWeightedCharge_;
        Float_t thirdleadingTrackPtWeightedCharge_;
        Float_t leadingMuonPtWeightedCharge_;
        Float_t leadingElectronPtWeightedCharge_;
        Float_t trackNumberWeightedJetPt_;
        Float_t chargeWeightedTrackId_;
        Float_t svChargeWeightedFlightDistance_;
        Float_t secondaryVertexCharge_;
        Float_t ipSignificanceLeadingTrack_;
        Float_t trueBJetId_;
        Float_t thereIsALeadingLepton_;
        Float_t thereIsALeadingMuon_;
        Float_t thereIsASecondaryVertex_;
    } jetChargeMvaStruct_;
    
};

#endif