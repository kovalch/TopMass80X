#ifndef Sample_h
#define Sample_h

#include <vector>

#include <TString.h>

#include "../../common/include/sampleHelpers.h"






/// Class defining a sample for processing keeping all information as needed
class Sample{
    
public:
    
    /// Specific type of sample as needed to be known for eg. plotting or Drell-Yan scale factor calculation
    /// Need to be also set for specific sample types after class definition
    enum SampleType{data, pseudodata, dyee, dymumu, dytautau, ttHbb, ttHother, ttbb, ttb, tt2b, ttcc, ttother, ttNoDilepton, ttZ, dummy};
    
    
    
    /// Default constructor
    Sample();
    
    /// Constructor for setting up a sample
    Sample(const TString& legendEntry,
           const int color,
           const double& crossSection,
           const double& crossSectionRelativeUp,
           const double& crossSectionRelativeDown,
           const std::vector<TString>& v_filename,
           const SampleType& sampleType =dummy);
    
    /// Default destructor
    ~Sample(){};
    
    
    /// Set sample legend entry for drawing
    void setLegendEntry(const TString& legendEntry);
    
    /// Return sample legend entry for drawing
    TString legendEntry()const;
    
    /// Set color of the sample
    void setColor(const int color);
    
    /// Return sample colour for drawing (needs to be identical for samples same legendEntry)
    int color()const;
    
    /// Return cross section corresponding to the sample
    double crossSection(const Systematic::Systematic& systematic)const;
    
    /// Set specific sample type
    void setSampleType(SampleType sampleType);
    
    /// Return the specific type of sample
    SampleType sampleType()const;
    
    /// Check if a specific filename is defined for this sample
    bool checkFilename(const TString& filename)const;
    
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
    
    /// Check whether sample contains the given root file
    bool containsFilenamesOfSample(const Sample& sample, const bool checkReweightedFiles =false)const;
    
    /// Get list of file names
    const std::vector<TString>& getFileNames()const;
    
    
    
private:
    
    /// Sample legend entry for drawing
    /// Samples will be ordered by legend entry and those with identical ones are merged in certain steps of further processing
    TString legendEntry_;
    
    /// Sample colour for drawing (needs to be identical for samples same legendEntry)
    int color_;
    
    /// Cross section corresponding to the sample
    double crossSection_;
    
    /// Cross section uncertainty up
    double crossSectionRelativeUp_;
    
    /// Cross section uncertainty down
    double crossSectionRelativeDown_;
    
    /// Specific type of sample as needed to be known for eg. plotting or Drell-Yan scale factor calculation
    SampleType sampleType_;
    
    /// Real final state of sample, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState_;
    
    /// Real systematic of sample, i.e. what should be used for given systematic (nominal or specific systematic)
    Systematic::Systematic systematic_;
    
    /// Path of the input root file
    TString inputFileName_;
    
    /// List of all file names (without channel prefix "channel_") for samples associated to this sample definition
    std::vector<TString> v_filename_;
};







#endif






