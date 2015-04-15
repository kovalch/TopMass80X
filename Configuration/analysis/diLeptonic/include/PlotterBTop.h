#ifndef PlotterBTop_h
#define PlotterBTop_h

#include <Rtypes.h>

class TCanvas;
class TLegend;
class TH1;
class TH2;
class TH2D;
class Output;
class TDirectory;

#include "TUnfoldBinning.h"
#include "SamplesFwd.h"

class PlotterBTop{
    
public:
    /// Constructor
    PlotterBTop(const Samples& samples, const double lumi, const double topxsec, const TString& specname);
    
    /// Destructor
    ~PlotterBTop(){};
    
    /// Set options specific for variables according to list in NameList
    void setOptions(const std::vector<TString>);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
     /// Prepare canvas and legend
     TLegend* setLegend();
    
private:
    void clearMemory();
    void clearMemoryPerSystematicChannel();
    
    ///Prepare Histograms
     void prepareHistograms(const std::vector<Sample>& v_sample);
     
     /// Write canvas
     void writeCanvas(TCanvas* cnavas, const TString name);
     
     const Samples& samples_;
     const double lumi_;
     const double topxsec_;
     const TString& specname_;
     
     TString plotName_;
     TString plotTitle_;
     TString plotUnits_;
     TString plotsFolder_;
     
    ///Control plots
    int cpNBins_;
    Double_t R1_;
    Double_t R2_;
    std::vector<TH1D* > v_SampleHist_; // v-sample , Hist 
    void writePlotCP(const std::vector<Sample>& v_sample,const std::vector<TH1D* >& v_SampleHist);
    TH2D* hRecoGen_;
    std::vector<Double_t> recogenBins_;
    
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
//     
    int genBin_(const float& val1,const float& val2);
    int recoBin_(const float& val1,const float& val2);
    void runUnfolding(const TH2* histMigration,const TH1* histInput,
                      const TH1* histBgr,const TH1* histBgrUo,
                      TH1*& histUf);
    
    ///Closure test
    TH1* histSR_;
    TH1* histUf_;
    TH1* histSG_;
//     TH1* histRWSG_;
//     double rewTopPtUp(double pt);
//     double rewTopPtDown(double pt);
//     double rewTopRapidityUp(double y);
//     double rewTopRapidityDown(double y);
//     void writePlotCT(TH1* histSG,TH1* histRW,TH1* histUf);
    
    /// Write EPS plots
    void writePlotEPS();
    TH1* histGen_;
    TH1* histPurity_;
    TH1* histStability_;
    TH1* histRecoGen_;
    
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

    std::vector<Double_t>  coarseBins_;
    std::vector<Double_t>  fineBins_;
    int uReco_;
    int oReco_;
    int uTrue_;
    int oTrue_;
    
    Double_t cut_;
    std::vector<Double_t>  cutCoarseBins_;
    std::vector<Double_t>  cutFineBins_;
  
    /// Tree branches
    Int_t entry_,entry0_;
    Float_t eventWeight_, trueLevelWeight_, trueLevelWeight0_;
    Int_t isTopGen_, isKinReco_;
    Float_t top_pt_;
    Float_t topbar_pt_;
    Float_t mlblbmet_;
    
    Float_t gen_top_pt_;
    Float_t gen_topbar_pt_;
    Float_t gen_mlblbmet_;
    
    Float_t gen0_top_pt_;
    Float_t gen0_topbar_pt_;
    Float_t gen0_mlblbmet_;
    
    Float_t branchVal_;
    Float_t branchValGen_;
    Float_t branchValGen0_;

    
    
    /// Histograms for Unfolding
    
//     TH1* histMCReco_;
//     TH2* histMCGenRec_;
    
    
    
    
};







#endif







