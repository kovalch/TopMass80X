#ifndef AnalyzerPlayground_h
#define AnalyzerPlayground_h

#include <vector>
#include <map>

class TString;
class TH1;

#include "AnalyzerBase.h"

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








/// Playground class, test here whatever you want to test
class AnalyzerPlayground : public AnalyzerBase{
    
public:
    
    /// Constructor
    AnalyzerPlayground(const std::vector<TString>& selectionStepsNoCategories,
                       const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                       const JetCategories* jetCategories =0);
    
    /// Destructor
    ~AnalyzerPlayground(){}
    
    
    
private:
    
    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
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







