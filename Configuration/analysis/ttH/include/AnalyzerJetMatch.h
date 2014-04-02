#ifndef AnalyzerJetMatch_h
#define AnalyzerJetMatch_h

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
};




#endif







