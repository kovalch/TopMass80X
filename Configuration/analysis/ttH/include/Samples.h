#ifndef Samples_h
#define Samples_h

#include <vector>
#include <map>

#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class RootFileReader;
class TString;



/// Class defining a sample for processing keeping all information as needed
class Sample{

public:

    /// Specific type of sample as needed to be known for eg. plotting or Drell-Yan scale factor calculation
    enum SampleType{data, dyee, dymumu, dytautau, ttHbb, ttHother, ttbb, ttb, ttother, dummy};



    /// Default constructor
    Sample();

    /// Constructor for setting up a sample
    Sample(const TString legendEntry, const int color, const double crossSection, const SampleType sampleType=dummy);

    /// Default destructor
    ~Sample(){};



    /// Return sample legend entry for drawing
    TString legendEntry()const;

    /// Return sample colour for drawing (needs to be identical for samples same legendEntry)
    int color()const;

    /// Return cross section corresponding to the sample
    double crossSection()const;

    /// Return the specific type of sample
    SampleType sampleType()const;



    /// Set real final state of sample, ie. only "ee", "emu", "mumu", but not "combined"
    void setFinalState(const Channel::Channel& channel);

    /// Get real final state of sample, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState()const;

    /// Set real systematic assigned to this sample, i.e. either nominal or specific systematic
    void setSystematic(const Systematic::Systematic& systematic);

    /// Get real systematic assigned to this sample, i.e. either nominal or specific systematic
    Systematic::Systematic systematic()const;

    /// Set the path of the input root file
    void setInputFile(const TString& inputFileName);

    /// Return the path of the input root file
    TString inputFile()const;



private:

    /// Sample legend entry for drawing
    /// Samples will be ordered by legend entry and those with identical ones are merged in certain steps of further processing
    TString legendEntry_;

    /// Sample colour for drawing (needs to be identical for samples same legendEntry)
    int color_;

    /// Cross section corresponding to the sample
    double crossSection_;

    /// Specific type of sample as needed to be known for eg. plotting or Drell-Yan scale factor calculation
    SampleType sampleType_;

    /// Real final state of sample, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState_;

    /// Real systematic of sample, i.e. what should be used for given systematic (nominal or specific systematic)
    Systematic::Systematic systematic_;

    /// Path of the input root file
    TString inputFileName_;
};





/// Class for administration of all samples, all dilepton channels and all systematics
class Samples{

public:

    /// Default constructor
    Samples();

    /// Constructor setting up samples
    Samples(const std::vector<Channel::Channel>& v_channel, const std::vector<Systematic::Systematic>& v_systematic);

    /// Default destructor
    ~Samples(){};



    /// Get map containing all samples per dilepton analysis channel and per systematic
    const SystematicChannelSamples& getSystematicChannelSamples()const;

    /// Get all samples of specific dilepton analysis channel and specific systematic
    const std::vector<Sample>& getSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)const;



private:

    /// Add samples for specific dilepton analysis channel and specific systematic
    void addSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic);

    /// Place where to define the samples as they will be used in the analysis
    std::vector<std::pair<TString, Sample> > setSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic);
    
    
    /// Combines samples with all corresponding input files
    /// Input files are identified by name patterns which need to be contained in the input file name
    std::vector<std::pair<TString, Sample> > samplesByNamePatterns(const std::vector<TString>& v_filename,
                                                                   const std::vector<std::pair<Sample, std::vector<TString> > >& v_sampleNamepatternsPair);
    
    /// Assign options to each sample via its filename
    std::vector<Sample> setSampleOptions(const Systematic::Systematic& systematic, const std::vector<std::pair<TString, Sample> >& v_filenameSamplePair);

    /// Order samples by their legend
    /// when a legend already exists, the sample is moved directly behind it
    void orderByLegend(std::vector<Sample>& v_sample);

    /// Assign the real final state to each sample, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel assignFinalState(const TString& filename);

    /// Assign the real systematic to each sample, i.e. what should be used for given systematic (nominal or specific systematic)
    /// and modify filename accordingly
    Systematic::Systematic assignSystematic(TString& filename, const Systematic::Systematic& systematic);



    /// Map containing all samples per dilepton analysis channel and per systematic
    SystematicChannelSamples m_systematicChannelSample_;
};





namespace Tools{

    /// Get the luminosity weight for a specific sample
    double luminosityWeight(const Sample& sample, const double luminosity, RootFileReader* fileReader);
}






#endif // Samples_h






