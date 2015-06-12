#ifndef AnalyzerMvaEventClassification_h
#define AnalyzerMvaEventClassification_h

#include <vector>
#include <map>

class TString;
class TH1;
class TH2D;

#include "AnalyzerBase.h"
#include "MvaVariablesTopJets.h"

class JetCategories;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinematicReconstructionSolutions;
class MvaVariablesTopJets;
class MvaReaderBase;
namespace tth{
    class RecoLevelWeights;
    class GenLevelWeights;
    class RecoObjectIndices;
    class GenObjectIndices;
}









/// Class for basic histograms that are filled simultaneously for any step
class AnalyzerMvaEventClassification : public AnalyzerBase{

public:

    /// Constructor
    AnalyzerMvaEventClassification(const std::vector<TString>& selectionStepsNoCategories,
                                   const std::vector<TString>& stepsForCategories = std::vector<TString>(),
                                   const JetCategories* jetCategories =0);

    /// Destructor
    ~AnalyzerMvaEventClassification(){}



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
    
    
    /// Map holding for each step (including category) a map of trainings and corresponding MvaReaders
    std::map<TString, std::map<TString, MvaReaderBase*> > m_m_mvaWeight_;
};




#endif







