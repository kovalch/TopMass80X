#ifndef HiggsAnalysis_h
#define HiggsAnalysis_h

#include <Rtypes.h>

class TTree;

#include "analysisHelpers.h"
#include "analysisStructsFwd.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"

class MvaTreeHandlerBase;
class AnalyzerBase;
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
    
    
    
    /// Set up all analysers of type AnalyzerBase
    void SetAllAnalyzers(std::vector<AnalyzerBase*> v_analyzer);
    
    /// Set up all tree handlers of type MvaTreeHandlerBase
    void SetAllTreeHandlers(std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler);
    
    
    
private:
    
    
    /// Select events from Higgs signal samples which need to be removed due to generator selection
    bool failsHiggsGeneratorSelection(const int higgsDecayMode)const;
    
    /// Select events from Top signal that satisfy flavour of the additional jets
    bool failsAdditionalJetFlavourSelection(const Long64_t& entry)const;
    
    
    
    /// Returns a vector of indices of gen jets which are in acceptance
    std::vector<int> genJetIndices(const VLV& allGenJets, const TopGenObjects& topGenObjects)const;
    
    /// Create vector of size of gen jets, and assign for each element a vector of indices of associated B hadrons
    std::vector<std::vector<int> > matchBhadronsToGenJets(const std::vector<int>& genJetIndices, const VLV& allGenJets, 
                                                          const TopGenObjects& topGenObjects)const;
    
    /// Create vector of size of gen jets, and assign for each element a vector of indices of associated C hadrons
    std::vector<std::vector<int> > matchChadronsToGenJets(const std::vector<int>& genJetIndices, const VLV& allGenJets, 
                                                          const TopGenObjects& topGenObjects)const;
    
    /// Returns vector of indices of gen jets containing B hadrons
    std::vector<int> genBjetIndices(const std::vector<std::vector<int> >& genJetBhadronIndices)const;
    
    /// Returns vector of indices of gen jets containing C hadrons, but no B hadrons
    std::vector<int> genCjetIndices(const std::vector<std::vector<int> >& genJetBhadronIndices,
                                    const std::vector<std::vector<int> >& genJetChadronIndices)const;
    
    /// Get indices of generated (anti-)b jet by its mother particle ID ("signed" mother PDG ID, "+" is b, "-" is anti-b)
    /// Requires that exactly one (anti-)b quark is found for the given ID
    /// Returns index of corresponding gen jet, or negative value for indicating cases no jet could be identified
    int genBjetIndex(const TopGenObjects& topGenObjects, const int pdgId)const;
    
    
    /// Match generated jet with given index to reco jets
    /// Returns the index of the matched reco jet, or negative value for indicating cases no jet could be identified
    int matchRecoToGenJet(const std::vector<int>& jetIndices, const VLV& jets, const int genJetIndex, const VLV& genJets)const;
    
    /// Create vector of size of gen jets,
    /// and assign for each element the index of the matched reco jet, or negative value for indicating cases no jet could be identified
    std::vector<int> matchRecoToGenJets(const std::vector<int>& jetIndices, const VLV& jets,
                                        const std::vector<int>& genJetIndices, const VLV& allGenJets)const;
    
    
    
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
    
    
    
    /// All analysers of type AnalyzerBase
    std::vector<AnalyzerBase*> v_analyzer_;
    
    /// All tree handlers of type MvaTreeHandlerBase
    std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler_;
};





#endif




