#ifndef EventYields_h
#define EventYields_h

#include <utility>
#include <vector>

class TString;

#include "../../common/include/sampleHelpers.h"

class RootFileReader;
class Sample;
class Samples;
class DyScaleFactors;





class EventYields{
    
public:
    
    /// Constructor for producing event yield tables
    EventYields(const char* outputDirectory, const Samples& samples);
    
    /// Default destructor
    ~EventYields(){};
    
    
    
private:
    
    /// Produce the yields
    void produceYields(const char* outputDirectory, const Samples& samples)const;
    
    /// Write the yields to txt files, either without or with additional corrections (e.g. Drell-Yan scaling)
    void writeYields(const char* outputDirectory,
                     const Samples& samples,
                     const std::pair<TString, TString>& nameStepPair,
                     const bool useCorrections =false)const;
    
    
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
};





#endif





