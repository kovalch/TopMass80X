#include "VariablesBase.h"





VariablesBase::VariablesBase():
eventWeight_(VariableFloat(name_eventWeight_))
{}



VariablesBase::VariablesBase(const double& eventWeight):
eventWeight_(VariableFloat(name_eventWeight_))
{
    // Event weight
    eventWeight_.value_ = eventWeight;
}



void VariablesBase::clearVariables(std::vector<VariablesBase*>& v_Variables){
    for(auto& variables : v_Variables) delete variables;
    v_Variables.clear();
}







