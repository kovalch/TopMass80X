#include "MvaVariablesBase.h"





MvaVariablesBase::MvaVariablesBase():
eventWeight_(MvaVariableFloat(name_eventWeight_))
{}



MvaVariablesBase::MvaVariablesBase(const double& eventWeight):
eventWeight_(MvaVariableFloat(name_eventWeight_))
{
    // Event weight
    eventWeight_.setValue(eventWeight);
}



void MvaVariablesBase::clearVariables(std::vector<MvaVariablesBase*>& v_mvaVariables){
    for(auto& mvaVariables : v_mvaVariables) delete mvaVariables;
    v_mvaVariables.clear();
}







