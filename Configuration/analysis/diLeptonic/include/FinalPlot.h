#ifndef FinalPlot_h
#define FinalPlot_h

#include <Rtypes.h>
#include <TCanvas.h>

class TCanvas;
//class TLegend;
//class TH1;
//class TH2;
//class TDirectory;
class Samples;

#include "SamplesFwd.h"

class FinalPlot{
    
public:
    
    /// Constructor
    FinalPlot(const std::vector<Channel::Channel>& v_channel,
              const std::vector<Systematic::Systematic>& v_systematic,
              const double lumi, const double topxsec);
    
    /// Destructor
    ~FinalPlot(){};
    
    /// Set options specific for variables according to list in NameList
    void setOptions(const std::vector<TString>);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots(const TString& fileType);
    
    /// Print core part of table with syst. and stat. errors for latex input.
    void printSystTableCore(const TString& name1, const TString& name2, const int& nBins,
                            const std::vector<std::vector<double> >& vv, const std::vector<std::vector<TString> >& vvLinks,
                            const TString& saveFile);
    
    
private:
    
    const std::vector<Channel::Channel>& v_channel_;
    const std::vector<Systematic::Systematic>& v_systematic_;
    const double lumi_;
    const double topxsec_;
    int nD_;
    
    std::vector<TString> v_plotName_;
    std::vector<TString> v_plotTitle_;
    std::vector<TString> v_plotUnits_;
    std::vector<TString> v_plotYunits_;
    int nBinsVar2_;
    
    void clearMemory();
    void clearMemoryPerSystematicChannel();
//     std::vector<TString> v_plotTitle_;
//     std::vector<TString> v_plotUnits_;
//     std::vector<std::vector<Double_t> >  v_coarseBins_;
    
   // ///Create canvas for plotting 
   // TCanvas* setCanvas();
    
    
    
};







#endif







