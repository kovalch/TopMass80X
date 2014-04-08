#ifndef GlobalScaleFactors_h
#define GlobalScaleFactors_h

#include <vector>
#include <map>

class TString;

#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class RootFileReader;
class DyScaleFactors;
class TtbbScaleFactors;




/// Class for handling global scale factors, i.e. which are valid per sample
class  GlobalScaleFactors{
    
public:
    
    struct ScaleFactorStruct{
        ScaleFactorStruct(const TString& step, const bool dyCorrection =false, const bool ttbbCorrection =false);
        ~ScaleFactorStruct(){}
        
        bool scaleFactorsExist(const TString& step, const bool dyCorrection, const bool ttbbCorrection)const;
        
        TString step_;
        bool dyCorrection_;
        bool ttbbCorrection_;
        
        bool anyCorrectionApplied_;
        
        SystematicChannelFactors systematicChannelScaleFactors_;
    };
    
    
    
public:
    
    GlobalScaleFactors(const std::vector<Channel::Channel>& v_channel,
                       const std::vector<Systematic::Systematic>& v_systematic,
                       const double& luminosityInInverseFb =1.,
                       const bool dyCorrection =false,
                       const bool ttbbCorrection =false);
    ~GlobalScaleFactors(){}
    
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step,
                                                           const bool dyCorrection, const bool ttbbCorrection)const;
    
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step)const;
    
    
    
private:
    
    
    
    DyScaleFactors* dyScaleFactors_;
    
    TtbbScaleFactors* ttbbScaleFactors_;
    
//    std::vector<ScaleFactorStruct> v_scaleFactorStruct_;
    
    SystematicChannelFactors m_luminosityWeight_;
};







#endif






