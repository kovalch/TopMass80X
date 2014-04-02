#ifndef AnalyzerDijet_h
#define AnalyzerDijet_h

#include <vector>
#include <string>
#include <map>

class TTree;
class TSelectorList;
class TH1;
class TString;

#include "JetCategories.h"
#include "AnalyzerBase.h"
#include "../../common/include/storeTemplate.h"
#include "../../common/include/classesFwd.h"

class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinRecoObjects;
namespace tth{
    class RecoObjectIndices;
    class GenObjectIndices;
}

class MvaTopJetsVariables;
class MvaReaderBase;



/// Class that analyzes all b-jet pairs from the input jet collection
/// In addition produces other plots about input jets, their origin, other properties
class AnalyzerDijet : public AnalyzerBase{

public:


    /// Constructor for given jet categories
    AnalyzerDijet(const char* mva2dWeightsFile, const std::string& corName, const std::string& swpName,
                  const std::vector<TString>& selectionStepsNoCategories,
                  const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                  const JetCategories* jetCategories =0, bool doHadronMatchingComparison = false);

    /// Empty destructor
    ~AnalyzerDijet(){};

    /// Find index of genJet corresponding to the specified reco jet. Returns -1 if not found
    int genJetIdOfRecoJet(const LV& recoJet, const VLV& genJets, const float dR_max=999.9);

    /// Get vector of indices of hadrons that are associted to the given gen jet
    std::vector<int> bHadIdsInGenJet(const int jetId, const std::vector<int>& hadJetIndices);

    /// Get vector of flavours of hadrons that are associted to the given gen jet
    std::vector<int> bHadFlavoursInGenJet(const int jetId, const std::vector<int>& hadJetIndices,
                                          const std::vector<int>& hadFlavours, const bool absFlavour = true);

    /// Whether index is in the vector of indices
    bool isInVector(const std::vector<int>& vector, const int id);
    bool isInVector(std::vector<std::pair<int,int> >& vector, const std::pair<int,int> id);

    /// Add to the vector only if it's not there already
    bool putUniquelyInVector(std::vector<int>& vector, const int id);
    bool putUniquelyInVector(std::vector<std::pair<int,int> >& vector, const std::pair<int,int> id);
    
    /// Creates a vector of all possible jet pairs from allJets except pairs containing any jet from jetsToIgnore
    std::vector<std::pair<int,int> > allPairsFromJets(const std::vector<int>& allJets, const std::vector<int>& jetsToIgnore);



private:

    /// Book histograms for one categoryId with given id and label
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book set of histograms for dijet mass
    void bookPairHistos(TH1* histo, std::map<TString, TH1*>& m_histogram, const TString& name);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const TString& step,
                            std::map<TString, TH1*>& m_histogram);
    
    /// Check how we additional b-jets in tt+bb events
    void checkAdditionalGenBJetAcceptance(const TopGenObjects& topGenObjects, const CommonGenObjects& commonGenObjects, 
                                          std::map<TString, TH1*>& m_histogram, const float jetPt_min, const float jetEta_max, const double weight);

    /// Analyze jet pairs of given jets for the given b-jets from top. Returns ration of correct pairs to wrong pairs
    float correctPairFraction(const VLV& allJets, const std::vector<int>& jetsId,
                              const std::vector<int>& bJetsId, const std::vector<double>& jetsBtagDiscriminant,
                              const std::vector<int>& topJetsId, const std::vector<int>& higgsJetsId,
                              const double weight, std::map<TString, TH1*>& m_histogram, std::string histoName, 
                              const bool fillAllCombinations = true, const double jetPt_threshold = 0.0, const int lowerHigher = 0, 
                              const std::vector<std::pair<int, int> > &pairsToIgnore = std::vector<std::pair<int, int> >(0) );
    
    /// Fill dijet mass for pairs of jets
    void fillDijetMassForPairs(const VLV& allJets, const std::vector<int>& jetsId, const std::vector<int>& higgsJetsId,
                               const std::vector<std::pair<int, int> > &jetPairs, const double weight, 
                               std::map<TString, TH1*>& m_histogram, std::string histoName, const bool normaliseWeight = false );
    
    /// Fill values for correct/wrong Higgs pairs separately
    void fillPairHistos(std::map<TString, TH1*>& m_histogram, const TString& name, const double value, const bool isCorrect, const double weight = 1.);

    /// Fill histograms about Gen/Reco matching: comparison of dR to true matching
    void fillGenRecoMatchingComparisonHistos(const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                             const VLV& bHadLVs, const std::vector<int>& bHadFlavour, const std::vector<int>& bHadJetIndex,
                                             const VLV& genJets, std::map<TString, TH1*>& m_histogram, const double weight);
    
    /// Returns list of jets pairs that are not from H according to the MVA
    std::vector<std::pair<int,int> > jetPairsFromMVA(std::map<TString, TH1*>& m_histogram, const tth::RecoObjectIndices& recoObjectIndices,
                                                     const tth::GenObjectIndices& genObjectIndices, const RecoObjects& recoObjects, 
                                                     const std::vector<int>& trueTopJetsId, const std::vector<int>& trueHiggsJetsId, 
                                                     const double weight);

    /// Checks whether two indices are among vector of pairs of indices
    bool areAmongPairs(const std::vector<std::pair<int, int> > &pairs, const int idx1, const int idx2);



    /// Struct holding the histograms for one jet category
    struct CatHistograms{
        /// Constructor
        CatHistograms(){};
        /// Destructor
        ~CatHistograms(){};

        /// The map with all the histograms for one jet category
        std::map<TString, TH1*> map_histo;
    };

    /// Histograms for all jet categories
    std::vector<CatHistograms> histograms_;

    /// MVA weights of correct dijet assignment for top system
    MvaReaderBase* weightsCorrect_;

    /// MVA weights of swapped dijet assignment for top system
    MvaReaderBase* weightsSwapped_;

    /// Whether to do comparison of dR and improved hadron-quark-jet matching
    bool doHadronMatchingComparison_;
    
    double sigJetPt_min_;
    double sigJetEta_max_;

};


#endif