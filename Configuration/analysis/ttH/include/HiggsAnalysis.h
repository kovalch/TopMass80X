#ifndef HiggsAnalysis_h
#define HiggsAnalysis_h

#include <Rtypes.h>

class TTree;

#include "analysisHelpers.h"
#include "analysisStructsFwd.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"

class MvaTreeHandlerBase;
class AnalyzerBaseClass;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinRecoObjects;
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
    
    
    
    /// Set up all analysers of type AnalyzerBaseClass
    void SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer);
    
    /// Set up all tree handlers of type MvaTreeHandlerBase
    void SetAllTreeHandlers(std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler);
    
    
    
private:
    
    
    /// Select events from Higgs signal samples which need to be removed due to generator selection
    bool failsHiggsGeneratorSelection(const int higgsDecayMode)const;
    
    /// Select events from Top signal that satisfy flavour of the additional jets
    bool failsAdditionalJetFlavourSelection(const Long64_t& entry)const;
    
    
    
    /// Get indices of generated b jet and anti-b jet for particles with given PDG ID
    /// Returns whether a unique solution is found, and only in this case indices are set unequal to -1
    bool getGenBjetIndices(int& genBjetIndex, int& genAntiBjetIndex,
                           const TopGenObjects& topGenObjects, const int pdgId);
    
    /// Match the two generated input jets to the reco jet collection
    bool matchRecoToGenJets(int& matchedBjetIndex, int& matchedAntiBjetIndex,
                            const std::vector<int>& jetIndices,
                            const VLV& jets,
                            const LV* genBjet, const LV* genAntiBjet);
    
    
    
    /// Return vector of pair of indices for dijet combinations, each pair ordered by jet charge
    tth::IndexPairs chargeOrderedJetPairIndices(const std::vector<int>& jetIndices,
                                                const std::vector<double>& jetCharges);

    
    
    /// Fill all analysers and histograms in one method
    void fillAll(const std::string& selectionStep,
                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                 const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                 const KinRecoObjects& kinRecoObjects,
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
    
    
    
    /// All analysers of type AnalyzerBaseClass
    std::vector<AnalyzerBaseClass*> v_analyzer_;
    
    /// All tree handlers of type MvaTreeHandlerBase
    std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler_;
};





#endif




