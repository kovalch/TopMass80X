#ifndef AnalyzerGenEvent_h
#define AnalyzerGenEvent_h

#include <map>
#include <vector>

class TString;
class TH1;

#include "AnalyzerBase.h"

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
class AnalyzerGenEvent : public AnalyzerBase{
    
public:
    
    /// Constructor
    AnalyzerGenEvent(const std::vector<TString>& selectionStepsNoCategories,
                     const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                     const JetCategories* jetCategories =0);
    
    /// Destructor
    ~AnalyzerGenEvent(){}
    
    
    
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
    
    
    
    /// Book histograms for the event-wise properties
    void bookEventHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book histograms for the hadron-genJet-recoJet matching of top and Higgs jets
    void bookMatchingHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book histos for b jets and their B hadron constituents
    void bookBhadronHistos(const TString& topOrHiggsName, const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Book histos for the two jets either from top system or from Higgs
    /// Only for events with both jets of the requested type separated
    void bookTopOrHiggsHistos(const TString& topOrHiggsName,
                              const TString& step,
                              std::map<TString, TH1*>& m_histogram);
    
    /// Book histos for all four jets from top system and Higgs
    /// Only for events with all four jets separated
    void bookTopAndHiggsHistos(const TString& step,
                               std::map<TString, TH1*>& m_histogram);
    
    
    
    /// Fill histograms for the event-wise properties
    void fillEventHistos(const CommonGenObjects& commonGenObjects, 
                            const tth::GenObjectIndices& genObjectIndices,
                            const double& weight, std::map<TString, TH1*>& m_histogram);
    
    /// Fill histograms for the hadron-genJet-recoJet matching of top and Higgs jets
    void fillMatchingHistos(const tth::GenObjectIndices& genObjectIndices,
                            const double& weight, std::map<TString, TH1*>& m_histogram);
    
    /// Fill histos for b jets and their B hadron constituents
    void fillBhadronHistos(const TString& topOrHiggsName,
                           const tth::GenObjectIndices& genObjectIndices,
                           const double& weight,
                           std::map<TString, TH1*>& m_histogram);
    
    /// Fill histos for the two jets either from top system or from Higgs
    /// Only for events with both jets of the requested type separated
    void fillTopOrHiggsHistos(const TString& topOrHiggsName,
                              const CommonGenObjects& commonGenObjects,
                              const tth::GenObjectIndices& genObjectIndices,
                              const double& weight,
                              std::map<TString, TH1*>& m_histogram);
    
    /// Fill histos for all four jets from top system and Higgs
    /// Only for events with all four jets separated
    void fillTopAndHiggsHistos(const CommonGenObjects& commonGenObjects,
                               const tth::GenObjectIndices& genObjectIndices,
                               const double& weight,
                               std::map<TString, TH1*>& m_histogram);
};




#endif







