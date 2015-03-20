#ifndef JetCharge_h
#define JetCharge_h

#include <vector>
#include <map>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class TString;
class TFile;
class TH1;

#include "AnalyzerBase.h"
#include "../../common/include/classesFwd.h"
#include "../../common/include/storeTemplate.h"

class JetCategories;
class MvaVariablesJetCharge;
class MvaReaderBase;
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


//FIXME: description
/// jet track reweighting class, test here whatever you want to test
class JetCharge{
    
public:
    
    /// Constructor
    JetCharge(const bool mvaCharge, const bool correction);
    
    /// Destructor
    ~JetCharge(){}
    
    // FIXME: description
    /// Momentum weighted jet charge calculation for a given squeezing parameter x
    double pWeightedCharge(const int jetIndex, const LV& recoJet,
                           const std::vector<int>& pfCandidateTrackIndex, const VLV& pfCandidates,
                           const std::vector<int>& pfCandidateCharge, const std::vector<int>& pfCandidateVertexId,
                           const double& x)const;
   
    /// MVA jet charge 
    double mvaCharge(const int jetIndex, const LV& jet, const RecoObjects& recoObjects);
    
    /// Vector of MVA jet charges for all jets
    std::vector<float> mvaCharges(const EventMetadata& eventMetadata,
                                  const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                  const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                  const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                                  const double& weight);
    
    /// Correct value of jet charge by applying a quantile mapping
    void quantileMappingCorrection(double& jetCharge)const;
    
    
    
private:
    
    /// Check if histogram with given name exists in file, and access clone in memory
    TH1* readHist(TFile* file, const TString& histname)const;
    
    // FIXME: Doxygen description of all variables
    
    MvaReaderBase* mvaReader_;
    
    const bool mvaCharge_;
    const bool correction_;
    
    TH1* histData_;
    TH1* histMc_;
    
    
    
    // FIXME: Remove and use generic ones
    TMVA::Reader* mvaReader2_;
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