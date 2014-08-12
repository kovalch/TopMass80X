#ifndef Correctors_h
#define Correctors_h

#include <string>

#include "classesFwd.h"

namespace ztop{
    class RecoilCorrector;
}






class MetRecoilCorrector{
    
public:
    
    /// Constructor
    MetRecoilCorrector(const std::string& fileData, const std::string& fileMC);
    
    /// Destructor
    ~MetRecoilCorrector();
    
    
    
    /// Apply correction to MET by modifying px and py
    double applyCorrection(LV& met, const LV& genZ, const LV& dilepton, const int nJet)const;
    
    
    
private:
    
    /// Pointer to the recoil corrector instance
    ztop::RecoilCorrector* const recoilCorrector_;
};






#endif





