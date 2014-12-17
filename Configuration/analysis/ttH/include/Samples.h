#ifndef Samples_h
#define Samples_h

#include <vector>
#include <map>
#include <utility>

class TString;

#include "Sample.h"
#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class GlobalScaleFactors;





/// Class for administration of all samples, all dilepton channels and all systematics
class Samples{
    
public:
    
    /// Default constructor
    Samples();
    
    /// Constructor setting up samples
    Samples(const TString& filelistDirectory,
            const std::vector<Channel::Channel>& v_channel,
            const std::vector<Systematic::Systematic>& v_systematic,
            const GlobalScaleFactors* globalScaleFactors =0);
    
    /// Default destructor
    ~Samples(){};
    
    
    
    /// Get map containing all samples per dilepton analysis channel and per systematic
    const SystematicChannelSamples& getSystematicChannelSamples()const;
    
    /// Get all samples of specific dilepton analysis channel and specific systematic
    const std::vector<Sample>& getSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)const;
    
    /// Set the class pointer to global scale factors to be associated to the Samples
    void setGlobalWeights(const GlobalScaleFactors* globalScaleFactors);
    
    /// Apply weights to samples based on luminosity and additional corrections if specified
    /// The selection step is extracted from the object name
    std::pair<SystematicChannelFactors, bool> globalWeights(const TString& objectname,
                                                            const bool dyCorrection,
                                                            const bool ttbbCorrection)const;
    
    /// Apply weights to samples based on all corrections which are set up for the GlobalScaleFactors
    /// The selection step is extracted from the object name
    std::pair<SystematicChannelFactors, bool> globalWeights(const TString& objectname)const;
    
    /// Return the used luminosity value in inverse pb, as it is stored in the GlobalScaleFactors
    double luminosityInInversePb()const;
    
    
    
    
private:
    
    /// Add samples for specific dilepton analysis channel and specific systematic
    void addSamples(const TString& filelistDirectory,
                    const Channel::Channel& channel,
                    const Systematic::Systematic& systematic);
    
    /// Set samples to be used at 8 TeV, and order them in the given order
    std::vector<std::pair<TString, Sample> > setSamples(const std::vector<TString>& v_filename,
                                                        const std::map<TString, Sample>& m_samples,
                                                        const std::vector<TString>& v_sampleIdentifier,
                                                        const bool hasPseudodata=false)const;
    
    
    
    /// Assign options to each sample via its filename
    /// For systematic variations of specific samples, set all other samples to nominal ones
    std::vector<Sample> setSampleOptions(const Systematic::Systematic& systematic,
                                         const std::vector<std::pair<TString, Sample> >& v_filenameSamplePair,
                                         const std::vector<std::pair<TString, Sample> >& v_filenameSamplePairNominal);
    
    /// Order samples by their legend
    /// when a legend already exists, the sample is moved directly behind it
    void orderByLegend(std::vector<Sample>& v_sample, const std::vector<TString>& v_legend)const;
    
    
    
    /// Map containing all samples per dilepton analysis channel and per systematic
    SystematicChannelSamples m_systematicChannelSample_;
    
    /// Pointer to the global scale factors
    const GlobalScaleFactors* globalScaleFactors_;    
};







#endif






