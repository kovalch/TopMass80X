#ifndef AnalyzerJetMatch_h
#define AnalyzerJetMatch_h

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
class AnalyzerJetMatch : public AnalyzerBase{
    
public:
    
    /// Constructor
    AnalyzerJetMatch(const std::vector<TString>& selectionStepsNoCategories,
                     const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                     const JetCategories* jetCategories =0);
    
    /// Destructor
    ~AnalyzerJetMatch(){}
    
    
    
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
    
    
    
    /// Book jet specific histos
    void bookJetHistos(const TString& whichSelection, const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book jet histos for all jets, and separately for top jets only and Higgs jets only
    void bookJetHistosInclExcl(const TString& whichSelection, const TString& whichJets,
                               const TString& step, std::map<TString, TH1*>& m_histogram);
    
    
    
    /// Fill jet specific histos
    void fillJetHistos(const TString& whichSelection,
                       const RecoObjects& recoObjects, const TopGenObjects& topGenObjects,
                       const int genIndex, const int recoIndex,
                       const tth::GenObjectIndices& genObjectIndices,
                       const double& weight,
                       std::map<TString, TH1*>& m_histogram);
    
    /// Fill jet histos for all jets, and separately for top jets only and Higgs jets only
    void fillJetHistosInclExcl(const TString& whichSelection, const TString& whichJets,
                               const LV& genJet, const LV& recoJet,
                               const double& weight,
                               std::map<TString, TH1*>& m_histogram);
    
    
    
    /// Check for the specified genJet is matched to the same reco jet as another one of the specified genJet collection
    bool isAmbiguous(const int genIndex, const std::vector<int>& genIndices, const std::vector<int>& closestRecoIndices)const;
};




#endif







