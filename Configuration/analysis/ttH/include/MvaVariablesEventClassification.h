#ifndef MvaVariablesEventClassification_h
#define MvaVariablesEventClassification_h


#include "MvaVariablesBase.h"

class RecoObjects;
namespace tth{
    class RecoObjectIndices;
}






class MvaVariablesEventClassification : public MvaVariablesBase{
    
public:
    
    /// Empty constructor
    MvaVariablesEventClassification();
    
    /// Constructor setting up input variables from physics objects
    MvaVariablesEventClassification(const double& btagDiscriminatorAverageTagged,
                                    const double& btagDiscriminatorAverageUntagged,
                                    const double& eventWeight);
    
    /// Destructor
    ~MvaVariablesEventClassification(){}
    
    /// Fill the MVA input structs for one event
    static MvaVariablesEventClassification fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                         const RecoObjects& recoObjects,
                                                         const double& eventWeight);
    
    
    
    // The variables needed for MVA
    
    // FIXME: describe each variable in doxygen
    /// Variables for MVA
    MvaVariableFloat btagDiscriminatorAverageTagged_;
    MvaVariableFloat btagDiscriminatorAverageUntagged_;
    
    
    
private:
    
    // The names associated to the variables
    
    static constexpr const char* name_btagDiscriminatorAverageTagged_ = "btagDiscriminatorAverageTagged";
    static constexpr const char* name_averageBtagDiscriminatorUntagged_ = "btagDiscriminatorAverageUntagged";
};






#endif







