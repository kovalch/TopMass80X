#ifndef GlobalScaleFactors_h
#define GlobalScaleFactors_h

#include <vector>

class TString;

#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class RootFileReader;
class LuminosityScaleFactors;
class DyScaleFactors;
class TtbbScaleFactors;





/// Class for handling global scale factors, i.e. which are valid per sample
class  GlobalScaleFactors{
    
public:
    
    /// Constructor setting up all requested scale factors for all specified channels and systematics
    GlobalScaleFactors(const std::vector<Channel::Channel>& v_channel,
                       const std::vector<Systematic::Systematic>& v_systematic,
                       const double& luminosityInInversePb =1.,
                       const bool dyCorrection =false,
                       const bool ttbbCorrection =false);
    
    /// Destructor
    ~GlobalScaleFactors(){}
    
    /// Access the scale factors for a given selection step by requesting explicitely which corrections to be used
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step,
                                                           const bool dyCorrection, const bool ttbbCorrection)const;
    
    /// Access the scale factors for a given selection step, using all corrections which are set up in constructor
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step)const;
    
    /// Return the used luminosity value in inverse pb
    double luminosityInInversePb()const;
    
    
    
private:
    
    /// The luminosity in inverse pb
    const double luminosityInInversePb_;
    
    /// Pointer to the luminosity scale factors
    LuminosityScaleFactors* luminosityScaleFactors_;
    
    /// Pointer to the Drell-Yan scale factors
    DyScaleFactors* dyScaleFactors_;
    
    /// Pointer to the tt+HF scale factors
    TtbbScaleFactors* ttbbScaleFactors_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* rootFileReader_;
};





#endif





