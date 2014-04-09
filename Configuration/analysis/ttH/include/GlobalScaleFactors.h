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
    
    /// Constructor setting up all requested scale factors for all specified channels and systematics
    GlobalScaleFactors(const std::vector<Channel::Channel>& v_channel,
                       const std::vector<Systematic::Systematic>& v_systematic,
                       const double& luminosityInInverseFb =1.,
                       const bool dyCorrection =false,
                       const bool ttbbCorrection =false);
    
    /// Destructor
    ~GlobalScaleFactors(){}
    
    /// Access the scale factors for a given selection step by requesting explicitely which corrections to be used
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step,
                                                           const bool dyCorrection, const bool ttbbCorrection)const;
    
    /// Access the scale factors for a given selection step, using all corrections which are set up in constructor
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step)const;
    
    
    
private:
    
    /// Pointer to the Drell-Yan scale factors
    DyScaleFactors* dyScaleFactors_;
    
    /// Pointer to the tt+HF scale factors
    TtbbScaleFactors* ttbbScaleFactors_;
    
    /// Map containing the luminosity weights for all samples
    SystematicChannelFactors m_luminosityWeight_;
};







#endif






