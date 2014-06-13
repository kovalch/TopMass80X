#ifndef SampleDefinitions_h
#define SampleDefinitions_h

#include <vector>
#include <map>

class TString;

class Sample;







/// Namespace for keeping sample definitions and sample ordering
namespace SampleDefinitions{
    
    /// All sample definitions for 8 TeV
    std::map<TString, Sample> samples8TeV();
    
    /// Select which samples for 8 TeV to use and in which order they should appear
    std::vector<TString> selectAndOrderSamples8TeV();
}







#endif






