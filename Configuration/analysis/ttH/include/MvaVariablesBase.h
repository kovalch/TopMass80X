#ifndef MvaVariablesBase_h
#define MvaVariablesBase_h

#include <string>

#include <Rtypes.h>

class TBranch;





// /// Templated struct which defines all relevant parameters for MVA input variables
// template<class T> struct MvaVariable{
//     MvaVariable(const char* name):
//         branch_(0), name_(name){}
//     
// public:
//     T value_;
//     TBranch* branch_;
//     
//     std::string name()const{return name_;}
//     char type()const{return type_;}
//     
// private:
//     const char* name_;
//     static constexpr char type_ = 'I';
// };



/// Struct which defines all relevant parameters for MVA input variables of type int
struct MvaVariableInt{
    MvaVariableInt(const char* name):
        value_(-999), valueFloat_(-999.F), branch_(0), name_(name){}
    
public:
    void setValue(const int value){value_ = value; valueFloat_ = static_cast<float>(value);}
    
    Int_t value_;
    Float_t valueFloat_; // Needed as newer root versions do not allow Int in TMVA reader anymore !?!
    TBranch* branch_;
    
    std::string name()const{return name_;}
    const char* type()const{return type_;}
    
private:
    const char* name_;
    static constexpr const char* type_ = "I";
};



/// Struct which defines all relevant parameters for MVA input variables of type float
struct MvaVariableFloat{
    MvaVariableFloat(const char* name):
        value_(-999.F), branch_(0), name_(name){}
    
public:
    void setValue(const float value){value_ = value;}
    
    Float_t value_;
    TBranch* branch_;
    
    std::string name()const{return name_;}
    const char* type()const{return type_;}
    
private:
    const char* name_;
    static constexpr const char* type_ = "F";
};



/// Base class for holding MVA input variables
class MvaVariablesBase{
    
public:
    
    /// Empty constructor
    MvaVariablesBase();
    
    /// Constructor setting up event weight
    MvaVariablesBase(const double& eventWeight);
    
    /// Destructor
    virtual ~MvaVariablesBase(){}
    
    
    
    /// Clear the MVA input variables, i.e. delete all pointers properly
    static void clearVariables(std::vector<MvaVariablesBase*>& v_mvaVariables);
    
    
    
    /// Event weight including all scale factors
    MvaVariableFloat eventWeight_;
    
    
    
private:
    
    static constexpr const char* name_eventWeight_ = "eventWeight";
};






#endif







