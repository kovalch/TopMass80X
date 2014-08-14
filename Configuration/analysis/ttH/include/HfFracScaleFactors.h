#ifndef hfFracScaleFactors_h
#define hfFracScaleFactors_h

#include <map>

class TString;

#include "../../common/include/sampleHelpers.h"
#include "Sample.h"

class RootFileReader;
class Samples;





class HfFracScaleFactors{

public:
    
    /// Constructor for producing Heavy-Flavour fraction scale factors
    HfFracScaleFactors(const Samples& samples, RootFileReader* const rootFileReader);
    
    /// Default destructor
    ~HfFracScaleFactors(){}
    
    
    
    /// Apply Heavy-Flavour fraction scale factor
    /// Returns 0 if it is not a sample that needs the scale factor and nothing needs to be done
    /// Returns -1 in case of no available scale factor for this step
    /// Returns +1 in case of successful application of scale factor
    int applyScaleFactor(double& weight,
                         const TString& fullStepname,
                         const Sample& sample,
                         const Systematic::Systematic& systematic)const;
    
    /// Normalising histogram
    void normalize ( TH1* histo )const;
    
private:
    
    /// Get Heavy-Flavour fraction scale factor for given selection step, systematic and channel
    const double& hfFracScaleFactor(const TString& step,
                                    const Systematic::Systematic& systematic,
                                    const Channel::Channel& channel,
                                    const Sample::SampleType& sampleType)const;
    
    
    /// Produce the Heavy-Flavour fraction scale factors
    void produceScaleFactors(const Samples& samples);
    
    /// Produce the Heavy-Flavour fraction scale factors for each selection step
    void produceScaleFactors(const TString& step, const Samples& samples);
    
    /// Struct to hold value-error pair
    struct ValErr {
        double val; double err;
        ValErr(double v, double e):val(v), err(e){}
        ValErr():val(1.), err(1.){}
            
    };
    
    /// Typedef for the map of scale factor value to the sample name
    typedef std::map<Sample::SampleType, ValErr> SampleTypeValueMap;
    
    /// Produce map of scale factors for each sample from fitting the histograms
    const std::vector<ValErr> getScaleFactorsFromHistos(const std::vector<TH1*> histos, const TString& step, 
                                                        const Channel::Channel channel)const;
    
    
    /// Typedef for the map containing the Heavy-Flavour fraction scale factors
    typedef std::map<TString, std::map<Systematic::Systematic, std::map<Channel::Channel, SampleTypeValueMap > > > HfFracScaleFactorMap;
    
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* const rootFileReader_;
    
    /// Map containing the Heavy-Flavour fraction scale factors
    HfFracScaleFactorMap m_hfFracScaleFactors_;
    
    /// Map containing the list of sample types with ids: to set specific combinations of samples to be fitted to data
    std::map<Sample::SampleType, int> sampleTypeIds_;
    
    /// Name base of the histogram used for the template fit
    TString histoTemplateName_;
};







#endif





