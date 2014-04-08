#ifndef dyScaleFactors_h
#define dyScaleFactors_h

#include <vector>
#include <map>

class TString;
class TH1;

#include "../../common/include/sampleHelpers.h"

class RootFileReader;
class Sample;
class Samples;





class DyScaleFactors{

public:
    
    /// Typedef for the map containing the Drell-Yan scale factors
    typedef std::map<TString, std::map<Systematic::Systematic, std::map<Channel::Channel, double> > > DyScaleFactorMap;
    
    
    
    /// Constructor for producing Drell-Yan scale factors
    DyScaleFactors(const Samples& samples);
    
    /// Default destructor
    ~DyScaleFactors(){}
    
    
    
    /// Apply Drell-Yan scale factor to histogram automatically depending on its name
    /// Returns 0 if it is not a Drell-Yan sample and nothing needs to be done
    /// Returns -1 in case of no available scale factor for this step
    /// Returns +1 in case of successful application of scale factor
    int applyScaleFactor(double& weight,
                         const TString& fullStepname,
                         const Sample& sample,
                         const Systematic::Systematic& systematic)const;
    
    
    
private:
    
    /// Get the map containing the Drell-Yan scale factors
    //const DyScaleFactorMap& getScaleFactors()const;
    
    /// Get Drell-Yan scale factor for given selection step, systematic and channel
    const double& dyScaleFactor(const TString& step,
                                const Systematic::Systematic& systematic,
                                const Channel::Channel& channel)const;
    
    
    
    /// Produce the Drell-Yan scale factors
    void produceScaleFactors(const Samples& samples);
    
    /// Produce the Drell-Yan scale factors for each selection step
    void produceScaleFactors(const TString& step, const Samples& samples);
    
    /// Print full information about all ingoing numbers (deactivated by default)
    void printFullInformation(const double dyScaleFactor_ee, const double dyScaleFactor_mumu, 
                              const double k_ee, const double k_mumu,
                              const double rOutIn_ee, const double rOutIn_mumu,
                              const double nIn_ee_data_loose, const double nIn_mumu_data_loose,
                              const double nIn_ee_data, const double nIn_mumu_data, const double nIn_emu_data,
                              const double nIn_ee_mc, const double nIn_mumu_mc,
                              const double nIn_ee_dy, const double nIn_mumu_dy,
                              const double nOut_ee_mc, const double nOut_mumu_mc,
                              const double nOut_ee_dy, const double nOut_mumu_dy,
                              const TString& step)const;
    
    
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    /// Map containing the Drell-Yan scale factors
    DyScaleFactorMap m_dyScaleFactors_;
};







#endif //dyScaleFactors_h





