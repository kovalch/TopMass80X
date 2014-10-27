#ifndef LuminosityScaleFactors_h
#define LuminosityScaleFactors_h

#include "SamplesFwd.h"

class RootFileReader;





class LuminosityScaleFactors{
    
public:
    
    /// Constructor setting up the luminosity scale factors for all samples
    LuminosityScaleFactors(const Samples& samples, const double& luminosityInInversePb, RootFileReader* const rootFileReader);
    
    /// Destructor
    ~LuminosityScaleFactors(){}
    
    /// Return the luminosity scale factor map
    const SystematicChannelFactors& scaleFactorMap()const;
    
    
    
private:
    
    /// Produce the luminosity scale factors
    void produceScaleFactors(const Samples& samples, const double& luminosityInInversePb);
    
    /// Get the luminosity weight, for a luminosity in inverse pb
    double luminosityWeight(const Sample& sample, const double& luminosityInInversePb)const;
    
    /// Calculate the luminosity weight
    double luminosityWeightPerInversePb(const Sample& sample)const;
    
    
    
    /// Map containing the luminosity weights for all samples
    SystematicChannelFactors m_weight_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* const rootFileReader_;
};






#endif






