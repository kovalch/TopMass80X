#ifndef SampleDefinitions_h
#define SampleDefinitions_h

#include <vector>
#include <set>
#include <map>

class TString;

class Sample;







/// Namespace for keeping sample definitions and sample ordering
namespace SampleDefinitions{
    typedef std::pair<int, std::set<TString> > ColorLegends;
    
    /// All sample definitions for 8 TeV
    std::map<TString, Sample> samples8TeV(const int mergeLevel);
    
    /// Select which samples for 8 TeV to use and in which order they should appear
    std::vector<TString> selectAndOrderSamples8TeV(const int pseudodata);
    
    /// All sample definitions for 13 TeV
    std::map<TString, Sample> samples13TeV(const int mergeLevel);
    
    /// Select which samples for 13 TeV to use and in which order they should appear
    std::vector<TString> selectAndOrderSamples13TeV(const int pseudodata);
    
    /// Vector of legend names (used to properly order samples)
    std::vector<TString> legendList(const std::map<TString, Sample>& samples, const std::vector<TString>& sampleIdentifiers);
    
    /// Whether pseudodata is used in the sample definitions
    bool usingPseudodata(const std::map<TString, Sample>& samples, const std::vector<TString>& sampleIdentifiers);
    
    /// List of merged legends for each broader legend-entry name
    std::map<TString, ColorLegends> mergedLegendEntries(const int mergeLevel = 1);
    
    /// Update legend entries of the merged samples
    void updateMergedLegendEntries(std::map<TString, Sample>& m_nameSample, const std::map<TString, ColorLegends>& m_legendLegends);
}







#endif






