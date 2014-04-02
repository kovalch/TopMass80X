#ifndef JetChargeAnalyzer_h
#define JetChargeAnalyzer_h

#include <map>
#include <vector>

class TString;
class TH1;

#include "AnalyzerBase.h"
#include "../../common/include/classesFwd.h"


class JetCategories;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinRecoObjects;
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
                      const JetCategories* jetCategories =0);

    /// Destructor
    ~AnalyzerJetCharge(){}

    /// Find index of genJet corresponding to the specified reco jet. Returns -1 if not found
    int genJetIdOfRecoJet(const LV& recoJet, const VLV& genJets, const float dR_max=999.9);
    
    /// Get vector of indices of hadrons that are associted to the given gen jet
    std::vector<int> bHadIdsInGenJet(const int jetId, const std::vector<int>& hadJetIndices);
    
    /// Get vector of flavours of hadrons that are associted to the given gen jet
    std::vector<int> bHadFlavoursInGenJet(const int jetId, const std::vector<int>& hadJetIndices,
					  const std::vector<int>& hadFlavours, const bool absFlavour = true);
	
    /// Whether index is in the vector of indices
    bool isInVector(const std::vector<int>& idVector, const int id);
    
    bool putUniquelyInVector(std::vector<int>& vector, const int id);
    
    struct MvaJetVariable
    {
        std::vector <float> longChargeJet_;
        std::vector <float> relChargeJet_;
        std::vector <int> leadingTrackCharge_;
        std::vector <float> leadingTrackPt_;
        std::vector <float> leadingTrackPtCharge_;
        std::vector <float> trueBJetPt_;
        std::vector <int> numTracks_;
        std::vector <int> trueBJetId_;
        std::vector<float> ptRatioTrackJet_;
        std::vector<bool> isMuonEvent_;
    } mvaStruct;
    
    TTree* mvaChargeTestTree;
    TTree* mvaChargeTrainTree;


private:

    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const TString& step,
                            std::map<TString, TH1*>& m_histogram);
};




#endif







