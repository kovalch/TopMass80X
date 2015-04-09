#ifndef Plotter_h
#define Plotter_h

#include <Rtypes.h>

class TCanvas;
class TLegend;
class TH1;
class TH2;
class Output;
class TDirectory;

#include "TUnfoldBinning.h"
#include "SamplesFwd.h"

class Plotter{
    
public:
    
    /// Constructor
    Plotter(const Samples& samples, const double lumi, double topxsec);
    
    /// Destructor
    ~Plotter(){};
    
    /// Set options specific for variables according to list in NameList
    void setOptions(const std::vector<TString>);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
private:
    
    void clearMemory();
    void clearMemoryPerSystematicChannel();
    
    ///Prepare Histograms
     void prepareHistograms(const std::vector<Sample>& v_sample);
     
     /// Write canvas
     void writeCanvas(TCanvas* cnavas, const TString name);
     
     const Samples& samples_;
     const double lumi_;
     double topxsec_;
     
     std::vector<TString> v_plotName_;
     std::vector<TString> v_plotTitle_;
     std::vector<TString> v_plotUnits_;
     
     TString plotsFolder_;
     
     
     
    ///Control plots
    std::vector<int> v_cpNBins_;
    std::vector<Double_t> v_R1_;
    std::vector<Double_t> v_R2_;
    std::vector<std::vector<TH1D* > > vv_SampleHist_; // v-nD_ , v-sample , Hist 
    std::vector<TH1* >  v_histRecoAllBins_;
    std::vector<std::vector<std::vector<TH1D* > > > vvv_SampleHist_; // v-nD_ , v-bin , v-sample , Hist    
    std::vector<std::vector<TH1D* > > vv_UnderflowSampleHist_; //v-nD_ , v-sample , Hist 
    std::vector<std::vector<TH1D* > > vv_OverflowSampleHist_;
    void writePlotCP(const std::vector<Sample>& v_sample,const std::vector<TH1D* >& v_SampleHist,const int& ind, const int binNum = -999);
    void writePlotCPAllBins(const std::vector<Sample>& v_sample ,const std::vector<TH1* >& v_SampleHist);
    
    ///Unfolding
    TUnfoldBinning* detectorBinning_;
    TUnfoldBinning* generatorBinning_;
    TUnfoldBinning* detectorDistribution_;
    TUnfoldBinning* generatorDistribution_;
    TH2* histMigration_;
    TH1* histBgrUo_;
    TH1* histBgr_;
    TH1* histData_;
    TH1* unfoldedData_;
    
    int genBin_(const std::vector<float>& val);
    int recoBin_(const std::vector<float>& val);
    void runUnfolding(const TH2* histMigration,const TH1* histInput,
                      const TH1* histBgr,const TH1* histBgrUo,
                      TH1*& histUf);
    
    ///Closure test
    TH1* histSR_;
    TH1* histUf_;
    TH1* histSG_;
    TH1* histRWSG_;
    double rewTopPtUp(double pt);
    double rewTopPtDown(double pt);
    double rewTopRapidityUp(double y);
    double rewTopRapidityDown(double y);
    void writePlotCT(TH1* histSG,TH1* histRW,TH1* histUf);
    
    TH1* histGen_;
    
    void writePlotEPSAllBins();
    TH1* histPurityAllBins_;
    TH1* histStabilityAllBins_;
    TH1* histRecoGenAllBins_;
    TH1* histEffAllBins_;
    TH1* histGenAllBins_;
    
    
    /// Write xSec plots
    void writePlotXSec(const TH1* hData,const TH1* hMC);
    std::vector<double> v_BR_;

    /// Access global scale factors
    const SystematicChannelFactors& scaleFactors(const TString& stepName);
    
    /// Map holding global scale factors for steps they were already accessed, to avoid re-calculating them
    std::map<TString, SystematicChannelFactors> m_stepFactors_;

    
    std::vector<std::vector<Double_t> >  v_coarseBins_;
    std::vector<std::vector<Double_t> >  v_fineBins_;
    std::vector<int> v_uReco_;
    std::vector<int> v_oReco_;
    std::vector<int> v_uTrue_;
    std::vector<int> v_oTrue_;
    
    /// Tree branches
    Int_t entry_,entry0_;
    Float_t eventWeight_, trueLevelWeight_, trueLevelWeight0_;
    Int_t isTopGen_, isKinReco_;
    std::vector<float> branchVals_;
    std::vector<float> branchValsGen_;
    std::vector<float> branchValsGen0_;
    
    /// Histograms for Unfolding
    
    TH1* histMCReco_;
    TH2* histMCGenRec_;
    
    
    
};







#endif







