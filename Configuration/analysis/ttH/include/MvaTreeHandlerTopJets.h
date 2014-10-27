#ifndef MvaTreeHandlerTopJets_h
#define MvaTreeHandlerTopJets_h

#include <vector>

class TString;
class TTree;

#include "MvaTreeHandlerBase.h"

class MvaVariablesBase;
class MvaTreePlotterBase;
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





/// Class for handling the trees of input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class MvaTreeHandlerTopJets : public MvaTreeHandlerBase{
    
public:
    
    /// Constructor for selection steps
    MvaTreeHandlerTopJets(const char* mvaInputDir,
                          const std::vector<TString>& selectionStepsNoCategories,
                          const std::vector<TString>& stepsForCategories =std::vector<TString>(),
                          const JetCategories* jetCategories =0);
    
    /// Destructor
    ~MvaTreeHandlerTopJets(){}
    
    
    
    /// Specify which type of plotter is needed for specific tree handler
    virtual MvaTreePlotterBase* setPlotter(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                           const bool separationPowerPlots =false)const;
    
    
    
private:
    
    /// Fill all variables for given selection step
    virtual void fillVariables(const EventMetadata& eventMetadata,
                               const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                               const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                               const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                               const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                               const double& weight, const TString& step,
                               std::vector<MvaVariablesBase*>& mvaVariables)const;
    
    /// Create and fill branches for TTree holding the input variables for MVA
    virtual void createAndFillBranches(TTree* tree, const std::vector<MvaVariablesBase*>& v_mvaVariables)const;
    
    /// Import all branches from TTree
    virtual void importBranches(TTree* tree, std::vector<MvaVariablesBase*>& mvaVariables)const;
};





#endif







