#ifndef JetChargeAnalyzer_h
#define JetChargeAnalyzer_h

#include <map>
#include <vector>

class TString;
class TH1;

#include "AnalyzerBase.h"
#include "../../common/include/classesFwd.h"
#include "../../common/include/storeTemplate.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


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



/// Class for basic histograms that are filled simultaneously for any step
class AnalyzerJetCharge : public AnalyzerBase{

public:
    
    /// Constructor
    AnalyzerJetCharge(const std::vector<TString>& selectionStepsNoCategories,
                      const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                      const JetCategories* jetCategories =0, const int debug=0);

    /// Destructor
    ~AnalyzerJetCharge(){}

    /// Check the decreasing order of a list (vector)
    bool checkDecreasingOrderOfAList(const VLV& vector);
    
    /// Check if the pfCandidate is matched to a selectedTrack and return the index of the selectedTrack
    int selTrackMatchToPfIndex(const LV& pfCandidate, const std::vector<LV>& selectedTracks, const int pfCharge, const std::vector<int>& selCharge, const int pfTrackJetIndex, const std::vector<int>& selectedTrackIndex, const std::vector <int>& ptOrderedSelTrackIdx);
    
    /// Pair jets with leptons in order to calculate the mlb. Returns a vector with three integers: first->sign of lepton (0 if mlb not succesfullly calculated), second->correctly matched jet (-1 if none), third->wrongly matched jet(-1 if none); successfull mlb means: one lepton+jet system under massThreshold, the other system above
    std::vector<int> leptonToJetMlbCalculator(const double massThreshold, const LV& lepton, const int leptonCharge, const VLV& recoJets, const std::vector<int>& jetIndices);
    
    /// Find index of genJet corresponding to the specified reco jet. Returns -1 if not found
    int genJetIdOfRecoJet(const LV& recoJet, const VLV& genJets, const float dR_max=999.9);
    
    /// Get vector of indices of hadrons that are associted to the given gen jet
    std::vector<int> bHadIdsInGenJet(const int jetId, const std::vector<int>& hadJetIndices);
    
    /// Get vector of flavours of hadrons that are associted to the given gen jet
    std::vector<int> bHadFlavoursInGenJet(const int jetId, const std::vector<int>& hadJetIndices,
					  const std::vector<int>& hadFlavours, const bool absFlavour = true);
	
    /// Whether index is in the vector of indices
    bool isInVector(const std::vector<int>& idVector, const int id);
    
    /// Check that an entry is only put once in the vector
    bool putUniquelyInVector(std::vector<int>& vector, const int id);
    
    /// Calculate the weight applying the track multiplicity in linear functional form (y = mx+n)
    double trackMultiplicityWeight(const double& m, const double& n, int jetIndex, const std::vector<int>& jetPfCandidateTrackIndex);
    
    /// Calculate the multiplicity per jet of a given collection
    unsigned int calculateMultiplicity(const std::vector<int>& collection, int jetIndex);
    
    /// p_{T} weighted jet charge calculation for a given squeezing parameter x (optimal value x = 0.8)
    double ptWeightedJetChargeX (const int jetId, const LV& recoJet, const double& x, const std::vector<int>& pfCandidateJetIndex, const VLV& pfCandidates, const std::vector<int>& pfCandidateCharge);
    
    /// Calculate product of all jet track multiplicity weights
    double trackMultiplicityWeightPerEvent (const std::vector<double>& jetWeight);
    
    struct MvaJetVariable
    {
        Float_t longChargeJet_;
        Float_t relChargeJet_;
        Float_t leadingTrackPtWeightedCharge_;
        Float_t leadingMuonPtWeightedCharge_;
        Float_t trackNumberWeightedJetPt_;
        Float_t chargeWeightedTrackId_;
        Float_t svChargeWeightedFlightDistance_;
        Float_t secondaryVertexCharge_;
        Float_t ipSignificanceLeadingTrack_;
        Float_t trueBJetId_;
        Float_t thereIsALeadingLepton_;
        Float_t thereIsALeadingMuon_;
        Float_t thereIsASecondaryVertex_;
    } mvaStruct_;
    
    TTree* mvaChargeTestTree_;
    TTree* mvaChargeTrainTree_;

    TMVA::Reader* reader_;
    TString case1_;

private:
    
    ///Debug boolean
    bool debug_;

    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book jet related (filled once per jet) histograms for given selection step
    void bookJetHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book pfCandidate related (filled once per pfCandidate) histograms for given selection step
    void bookPfCandidateHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book selected track related (filled once per selected track) histograms for given selection step
    void bookSelectedTrackHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book mva related (filled once per jet) histograms for given selection step
    void bookMvaHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book test and "other" related histograms for given selection step
    void bookOtherHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const EventMetadata& eventMetadata,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                            const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const TString& step,
                            std::map<TString, TH1*>& m_histogram);
};




#endif







