#ifndef HiggsAnalysis_h
#define HiggsAnalysis_h

#include <Rtypes.h>
#include <TString.h>

class TTree;

#include "analysisHelpers.h"
#include "analysisStructsFwd.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"

class MvaTreeHandlerBase;
class AnalyzerBase;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinematicReconstructionSolutions;
namespace tth{
    class GenLevelWeights;
    class RecoLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}




/// Class for the analysis step of the Higgs analysis, reading the nTuples
class HiggsAnalysis : public AnalysisBase{
    
public:
    
    /// Constructor
    HiggsAnalysis(TTree* =0);
    
    /// Destructor
    virtual ~HiggsAnalysis();
    
    /// Methods inherited from AnalysisBase
    virtual void Begin(TTree*);
    virtual void SlaveBegin(TTree*);
    virtual Bool_t Process(Long64_t entry);
    virtual void SlaveTerminate();
    virtual void Terminate();
    
    /// Class definition
    ClassDef(HiggsAnalysis, 0);
    
    
    
    /// ID for separating Higgs sample inclusive in Higgs decay via their decay mode
    void SetInclusiveHiggsDecayMode(const int inclusiveHiggsDecayMode);
    
    /// ID for separating ttbar samples via flavour of additional jets
    void SetAdditionalBjetMode(const int additionalBjetMode);
    
    /// Name of the reweighting
    void SetReweightingName(const TString& reweightingName);
    
    /// Slope of the reweighting
    void SetReweightingSlope(const double& reweightingSlope);
    
    /// Define for which processes to access genObjects already at beginning
    void SetGenStudies(const bool ttbb, const bool tth);
    
    
    
    /// Set up all analysers of type AnalyzerBase
    void SetAllAnalyzers(std::vector<AnalyzerBase*> v_analyzer);
    
    /// Set up all tree handlers of type MvaTreeHandlerBase
    void SetAllTreeHandlers(std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler);
    
    
    
private:
    
    
    /// Select events from Higgs signal samples which need to be removed due to generator selection
    bool failsHiggsGeneratorSelection(const int higgsDecayMode)const;
    
    /// Select events from Top signal that satisfy flavour of the additional jets
    bool failsAdditionalJetFlavourSelection(const int topDecayMode, const int additionalJetFlavourId)const;
    
    
    
    /// Return vector of pair of indices for dijet combinations, each pair ordered by jet charge
    tth::IndexPairs chargeOrderedJetPairIndices(const std::vector<int>& jetIndices,
                                                const std::vector<double>& jetCharges)const;
    
    
    /// Additional artifical weight to be applied for reweighting
    double reweightingWeight(const TopGenObjects& topGenObjects, const tth::GenObjectIndices& genObjectIndices)const;
    
    
    /// Select all reco object indices fulfilling selections
    void recoObjectSelection(std::vector<int>& allLeptonIndices,
                             std::vector<int>& leptonIndices, std::vector<int>& antiLeptonIndices,
                             int& leptonIndex, int& antiLeptonIndex,
                             int& leadingLeptonIndex, int& nLeadingLeptonIndex,
                             int& leptonXIndex, int& leptonYIndex,
                             std::vector<int>& jetIndices, tth::IndexPairs jetIndexPairs,
                             std::vector<int>& bjetIndices,
                             const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                             const Long64_t& entry)const;
    
    /// Select all gen object indices fulfilling selections
    void genObjectSelection(std::vector<int>& genJetIndices,
                            std::vector<std::vector<int> >& genJetBhadronIndices,
                            std::vector<int>& allGenBjetIndices, std::vector<int>& genBjetIndices,
                            std::vector<std::vector<int> >& genJetChadronIndices,
                            std::vector<int>& allGenCjetIndices, std::vector<int>& genCjetIndices,
                            int& genBjetFromTopIndex, int& genAntiBjetFromTopIndex,
                            int& genBjetFromHiggsIndex, int& genAntiBjetFromHiggsIndex,
                            const int higgsDecayMode, const int additionalJetFlavourId, const std::vector<int>& v_zDecayMode,
                            const TopGenObjects& topGenObjects)const;
    
    /// Match reco object indices to gen objects
    void matchRecoToGenObjects(std::vector<int>& genJetMatchedRecoBjetIndices,
                               std::vector<int>& genJetMatchedRecoCjetIndices,
                               int& matchedBjetFromTopIndex, int& matchedAntiBjetFromTopIndex,
                               int& matchedBjetFromHiggsIndex, int& matchedAntiBjetFromHiggsIndex,
                               const tth::GenObjectIndices& noRecoMatchGenObjectIndices,
                               const std::vector<int>& jetIndices, const VLV& jets,
                               const TopGenObjects& topGenObjects)const;
    
    
    /// Fill all analysers and histograms in one method
    void fillAll(const std::string& selectionStep,
                 const EventMetadata& eventMetadata,
                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                 const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                 const tth::GenObjectIndices& genObjectIndices, const tth::RecoObjectIndices& recoObjectIndices,
                 const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                 const double& defaultWeight)const;
    
    /// Book all histograms of all analysers for all steps in one method
    void bookAll();
    
    /// Clear all analysers in one method
    void clearAll();
    
    
    
    /// For a Higgs sample inclusive in decay, select H->bbbar or H->other decay (no separation for default value -999)
    int inclusiveHiggsDecayMode_;
    
    /// For a ttbar sample, select tt+bb, tt+b, or tt+other events (no separation for default value -999)
    int additionalBjetMode_;

    /// Name of the reweighting function
    TString reweightingName_;
    
    /// Slope of the reweighted shape
    double reweightingSlope_;
    
    /// Whether to access genObjects at first step for tt+bb
    bool genStudiesTtbb_;
    
    /// Whether to access genObjects at first step for ttH
    bool genStudiesTth_;
    
    
    
    /// All analysers of type AnalyzerBase
    std::vector<AnalyzerBase*> v_analyzer_;
    
    /// All tree handlers of type MvaTreeHandlerBase
    std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler_;
};





#endif




