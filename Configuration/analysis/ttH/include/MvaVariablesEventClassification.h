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
    MvaVariablesEventClassification(const int multiplicity_jets,
                                    const double& btagDiscriminatorAverage_tagged, const double& btagDiscriminatorAverage_untagged,
                                    const double& minDeltaR_jet_jet, const double& ptSum_jets_leptons,
                                    const int multiplicity_higgsLikeDijet15,
                                    const double& mass_higgsLikeDijet, const double& mass_higgsLikeDijet2,
                                    const double& eventWeight);
    
    /// Destructor
    ~MvaVariablesEventClassification(){}
    
    /// Fill the MVA input structs for one event
    static MvaVariablesEventClassification* fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                          const RecoObjects& recoObjects,
                                                          const double& eventWeight);
    
    
    
    // The variables needed for MVA
    
    // FIXME: describe each variable in doxygen
    /// Variables for MVA
    MvaVariableInt multiplicity_jets_;
    MvaVariableFloat btagDiscriminatorAverage_tagged_;
    MvaVariableFloat btagDiscriminatorAverage_untagged_;
    MvaVariableFloat minDeltaR_jet_jet_;
    MvaVariableFloat ptSum_jets_leptons_;
    MvaVariableInt multiplicity_higgsLikeDijet15_;
    MvaVariableFloat mass_higgsLikeDijet_;
    MvaVariableFloat mass_higgsLikeDijet2_;
    
    
    
private:
    
    // The names associated to the variables
    
    static constexpr const char* name_multiplicity_jets_ = "multiplicity_jets";
    static constexpr const char* name_btagDiscriminatorAverage_tagged_ = "btagDiscriminatorAverage_tagged";
    static constexpr const char* name_averageBtagDiscriminatorUntagged_ = "btagDiscriminatorAverage_untagged";
    static constexpr const char* name_minDeltaR_jet_jet_ = "minDeltaR_jet_jet";
    static constexpr const char* name_ptSum_jets_leptons_ = "ptSum_jets_leptons";
    static constexpr const char* name_multiplicity_higgsLikeDijet15_ = "multiplicity_higgsLikeDijet15";
    static constexpr const char* name_mass_higgsLikeDijet_ = "mass_higgsLikeDijet";
    static constexpr const char* name_mass_higgsLikeDijet2_ = "mass_higgsLikeDijet2";
};






#endif







