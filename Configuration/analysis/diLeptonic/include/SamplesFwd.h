#ifndef SamplesFwd_h
#define SamplesFwd_h

#include <vector>
#include <map>
#include <utility>

class TH1;
class TTree;

#include "../../common/include/sampleHelpers.h"





// Here only forward declaration to define typedefs
class Sample;

/// Storage type of all samples to be used in current analysis
/// These are all samples per dilepton analysis channel and per systematic      
typedef std::map<Systematic::Systematic, std::map<Channel::Channel, std::vector<Sample> > > SystematicChannelSamples;

/// Storage type for factors in association to samples, e.g. used for scale factors
typedef std::map<Systematic::Systematic, std::map <Channel::Channel, std::vector<double> > > SystematicChannelFactors;

/// Type for associating a specific histogram of the input file to the sample
typedef std::pair<Sample, TH1*> SampleHistPair;

/// Type for associating a tree of the input file to the sample
typedef std::pair<Sample, TTree*> SampleTreePair;





// Here only forward declaration to define typedefs
class Samples;





#endif






