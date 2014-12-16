#ifndef FinalPlot_h
#define FinalPlot_h

#include <Rtypes.h>

//class TCanvas;
//class TLegend;
//class TH1;
//class TH2;
//class TDirectory;
class Samples;

#include "SamplesFwd.h"

class FinalPlot{
    
public:
    
    /// Constructor
    FinalPlot(const Samples& samples, const double lumi, const double topxsec);
    
    /// Destructor
    ~FinalPlot(){};
    
    /// Set options specific for variables according to list in NameList
    void setOptions(const std::vector<TString>);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    const Samples& samples_;
    const double lumi_;
    const double topxsec_;
    int nD_;
    
    std::vector<TString> v_plotName_;
    TString plotsFolder_;
    
    void clearMemory();
    void clearMemoryPerSystematicChannel();
//     std::vector<TString> v_plotTitle_;
//     std::vector<TString> v_plotUnits_;
//     
//     std::vector<std::vector<Double_t> >  v_coarseBins_;
    
};







#endif







