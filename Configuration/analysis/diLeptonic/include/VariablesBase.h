#ifndef VariablesBase_h
#define VariablesBase_h

#include <string>

#include <Rtypes.h>

class TBranch;



/// Struct which defines all relevant parameters for TUnfold input variables of type int
struct VariableInt{
    VariableInt(const char* name):
        value_(-999), branch_(0), name_(name){}
    
public:
    Int_t value_;
    TBranch* branch_;
    
    std::string name()const{return name_;}
    const char* type()const{return type_;}
    
private:
    const char* name_;
    static constexpr const char* type_ = "I";
};



/// Struct which defines all relevant parameters for TUnfold input variables of type float
struct VariableFloat{
    VariableFloat(const char* name):
        value_(-999.F), branch_(0), name_(name){}
    
public:
    Float_t value_;
    TBranch* branch_;
    
    std::string name()const{return name_;}
    const char* type()const{return type_;}
    
private:
    const char* name_;
    static constexpr const char* type_ = "F";
};



/// Base class for holding TUnfold input variables
class VariablesBase{
    
public:
    
    /// Empty constructor
    VariablesBase();
    
    /// Constructor setting up event weight
    VariablesBase(const double& eventWeight);
    
    /// Destructor
    virtual ~VariablesBase(){}
    
    
    
    /// Clear the TUnfold input variables, i.e. delete all pointers properly
    static void clearVariables(std::vector<VariablesBase*>& v_Variables);
    
    
    
    /// Event weight including all scale factors
    VariableFloat eventWeight_;
    
    
    
private:
    
    static constexpr const char* name_eventWeight_ = "eventWeight";
};






#endif







