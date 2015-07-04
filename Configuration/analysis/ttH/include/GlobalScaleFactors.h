#ifndef GlobalScaleFactors_h
#define GlobalScaleFactors_h

#include <vector>

class TString;

#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class AnalysisConfig;
class RootFileReader;
class LuminosityScaleFactors;
class DyScaleFactors;
class HfFracScaleFactors;





/// Class for handling global scale factors, i.e. which are valid per sample
class  GlobalScaleFactors{
    
public:
    
    /// Constructor setting up all requested scale factors for all specified channels and systematics
    GlobalScaleFactors(const AnalysisConfig& analysisConfig,
                       const std::vector<Channel::Channel>& v_channel,
                       const std::vector<Systematic::Systematic>& v_systematic,
                       const bool dyCorrection =false,
                       const bool hfFracCorrection =false,
                       const bool writeToFile =false);
    
    /// Destructor
    ~GlobalScaleFactors(){}
    
    /// Access the scale factors for a given selection step by requesting explicitely which corrections to be used
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step,
                                                           const bool dyCorrection, const bool hfFracCorrection)const;
    
    /// Access the scale factors for a given selection step, using all corrections which are set up in constructor
    std::pair<SystematicChannelFactors, bool> scaleFactors(const Samples& samples, const TString& step)const;
    
    /// Return the used luminosity value in inverse pb
    double luminosityInInversePb()const;
    
    /// Return the relative uncertainty of the luminosity
    double luminosityUncertaintyRelative()const;
    
    
private:
    
    /// The luminosity in inverse pb
    const double luminosityInInversePb_;
    
    /// The relative luminosity uncertainty
    const double luminosityUncertaintyRelative_;
    
    /// Pointer to the luminosity scale factors
    LuminosityScaleFactors* luminosityScaleFactors_;
    
    /// Pointer to the Drell-Yan scale factors
    DyScaleFactors* dyScaleFactors_;
    
    /// Pointer to the tt+HF scale factors
    HfFracScaleFactors* hfFracScaleFactors_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* rootFileReader_;
};





#endif





