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
    std::pair <double, int> pWeightedCharge(const int jetIndex, const LV& recoJet,
                           const std::vector<int>& pfCandidateTrackIndex, const VLV& pfCandidates,
                           const std::vector<int>& pfCandidateCharge, const std::vector<int>& pfCandidateVertexId,
                           const double& x)const;
   
    /// MVA jet charge 
    double mvaCharge(const int jetIndex, const LV& recoJet, const RecoObjects& recoObjects)const;
    
    /// Vector of MVA jet charges for all jets
    std::vector<float> mvaCharges(const EventMetadata& eventMetadata,
                                  const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                  const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                  const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                                  const double& weight);
    
    /// Correct value of jet charge by applying a quantile mapping
    double quantileMappingCorrection(double& jetCharge)const;
    
    /// Return the correct charge value
    double jetChargeValue(const int jetIndex, const LV& recoJet,
                          const std::vector<int>& pfCandidateTrackIndex, const VLV& pfCandidates,
                          const std::vector<int>& pfCandidateCharge, const std::vector<int>& pfCandidateVertexId,
                          const double& x, const bool& isMc, const RecoObjects& recoObjects)const;
    
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
        mutable Float_t longChargeJet_;
        mutable Float_t relChargeJet_;
        mutable Float_t leadingTrackPtWeightedCharge_;
        mutable Float_t subleadingTrackPtWeightedCharge_;
        mutable Float_t thirdleadingTrackPtWeightedCharge_;
        mutable Float_t leadingMuonPtWeightedCharge_;
        mutable Float_t leadingElectronPtWeightedCharge_;
        mutable Float_t trackNumberWeightedJetPt_;
        mutable Float_t chargeWeightedTrackId_;
        mutable Float_t svChargeWeightedFlightDistance_;
        mutable Float_t secondaryVertexCharge_;
        mutable Float_t ipSignificanceLeadingTrack_;
        mutable Float_t trueBJetId_;
        mutable Float_t thereIsALeadingLepton_;
        mutable Float_t thereIsALeadingMuon_;
        mutable Float_t thereIsASecondaryVertex_;
    } jetChargeMvaStruct_;
    
};

#endif