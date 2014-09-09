#ifndef PlotterDiffXS_h
#define PlotterDiffXS_h

#include <vector>
#include <map>

#include <TString.h>

class TLegend;
class RootFileReader;
class TH1;
class TH2;
class TGraph;
class TGraphErrors;
class TPaveText;

#include "plotterHelpers.h"
#include "SamplesFwd.h"
#include "Sample.h"
#include "../../common/include/sampleHelpers.h"





class PlotterDiffXS{
    
public:
    
    /// Constructor
    PlotterDiffXS(const char* outputDir,
                  const Samples& samples,
                  const double luminosity);
    
    /// Destructor
    ~PlotterDiffXS(){};
    
    /// Set options specific for histogram under consideration
    void setOptions(const TString&, const TString&, const TString&, const TString&,
                    const int, const bool, const bool, const bool,
                    const double&, const double&, const double&, const double&);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    /// Pair of a legend entry and the histogram for the corresponding sample
    typedef std::pair<TString, TH1*> LegendHistPair;
    
    /// Access histogram under consideration from all samples, do global scalings
    bool prepareDataset(const std::vector<Sample>& v_sample,
                        const std::vector<double>& v_weight,
                        std::vector<SampleHistPair>& v_sampleHistPair,
                        const TString name);
    
    /// Calculate the response matrix
    void writeResponseMatrix(const Channel::Channel& channel, const Systematic::Systematic& systematic);

    /// Calculate the differential cross section
    void writeDiffXS(const Channel::Channel& channel, const Systematic::Systematic& systematic);
    
    /// Access global scale factors
    const SystematicChannelFactors scaleFactors(const bool dyScale = true, const bool hfScale = true);
    
    
    
    /// Set the style of the plot
    void setStyle(SampleHistPair& sampleHistPair);
    
    /// Set the style of the graph
    void setGraphStyle( TGraph* graph, Style_t marker = 21, Color_t markerColor = 1, Size_t markerSize = 1, 
                        Style_t line = 0, Color_t lineColor = 0, Size_t lineWidth = 1)const;

    /// Update histogram axis
    void updateHistoAxis(TH1* histo)const;
    
    /// Draw label for decay channel in upper left corner of plot
    void drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize =0.04)const;

    /// Draw official labels (CMS [Preliminary], luminosity and CM energy) above plot
    void drawCmsLabels(const int cmsprelim =1, const double& energy =8, const double& textSize =0.04)const;
    
    /// Draw purity/stability plots for the response matrices
    void drawPurityStability(TH2* histo2d, TString name)const;
    
    /// Draw signal significance label over the plot
    TPaveText* drawSignificance(TH1* signal, TH1* bkg, float Xmin,  float Xmax, float yOffset = 0.f, std::string sLabel ="", const int type=0)const;
    
    /// Calculates purity (type=0) and stability (type=1) curves for a 2D distribution
    TGraphErrors* purityStabilityGraph(TH2* h2d, const int type)const;
    
    /// Calculates binomial uncertainty of the subset/set ratio
    double uncertaintyBinomial(const double pass, const double all)const;

    /// Creates a vector of legend-histo pairs of specific type (0-data, 1-signal, 2-background)
    std::vector<LegendHistPair> legendHistPairsForSamples(const std::vector<Sample::SampleType> allowedSampleTypes, std::vector<SampleHistPair> samples)const;
    
    /// Sums all histograms in the stack
    TH1* sumOfHistos(const std::vector<LegendHistPair> histos, const TString name = "")const;
    
    /// Calculate differential cross section
    std::map<TString, TH1*> calculateDiffXS(const TH1* h_data, const TH1* h_signal, const TH1* h_bkg, const TH1* h_signal_noWeight, 
                                            const TH1* h_signal_gen, const TH1* h_ttbar_gen)const;
    
    
    

    /// Output folder name
    const char* outputDir_;
    
    /// Samples to be analysed
    const Samples& samples_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    
    
    /// Sample-Histogram pairs for reco level quantity to plot
    std::vector<SampleHistPair> v_sampleHistPair_;
    std::vector<SampleHistPair> v_sampleHistPair_noWeight_;
    
    /// Sample-Histogram pairs for gen level quantity to plot
    std::vector<SampleHistPair> v_sampleHistPairGen_;
    
    /// Sample-Histogram pairs for gen level events
    std::vector<SampleHistPair> v_sampleHistPairGenEventBased_;
    
    
    /// Map holding global scale factors for steps they were already accessed, to avoid re-calculating them
    std::map<TString, SystematicChannelFactors> m_stepFactors_;
    
    
    
    /// Name of histogram under consideration
    TString name_;
    TString nameGen_;
    TString nameGenEventBased_;
    
    /// Options for the histogram under consideration
    int mode_, rebin_;
//     bool stackToNEntries_;
    double rangemin_, rangemax_, ymin_, ymax_;
//     std::vector<double> XAxisbins_, XAxisbinCenters_;
    std::vector<Sample::SampleType> sampleTypesData_;
    std::vector<Sample::SampleType> sampleTypesSignal_;
    std::vector<Sample::SampleType> sampleTypesBackground_;
    std::vector<Sample::SampleType> sampleTypesTtbar_;
    TString YAxis_, XAxis_;
    bool logX_, logY_; // The variable logX_ is not used at all...
    
    /// Data for cross section calculation
    double luminosity_;
    struct ValueError {
        double v; double e;
        ValueError(double val, double err):v(val), e(err){}
        ValueError():v(1.), e(1.){}
    };
};







#endif







