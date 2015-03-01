#ifndef TreeHandlerTTBar_h
#define TreeHandlerTTBar_h

#include <vector>

class TString;
class TTree;

#include "TreeHandlerBase.h"

class VariablesBase;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinematicReconstructionSolutions;
namespace ttbar{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}





/// Class for handling the trees of input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class TreeHandlerTTBar : public TreeHandlerBase{
    
public:
    
    /// Constructor for selection steps
    TreeHandlerTTBar(const char* inputDir,
                          const std::vector<TString>& selectionStepsNoCategories);
    
    /// Destructor
    ~TreeHandlerTTBar(){}
    
    
    
    
private:
    
    /// Fill all variables for given selection step
    virtual void fillVariables(const EventMetadata& eventMetadata, 
                               const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects,
                               const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                               const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                               const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                               const double& weight, const TString& step,
                               std::vector<VariablesBase*>& variables)const;

    /// Book branches for TTree variables (dummy method, override in inherited TreeHandler)
    virtual void bookBranches(TTree* tree, VariablesBase* const variables_)const;
    
    /// Fill branches for TTree variables
    virtual void fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables);
    
    /// Import all branches from TTree
    virtual void importBranches(TTree* tree, std::vector<VariablesBase*>& variables)const;
};





#endif







