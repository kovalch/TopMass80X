#ifndef Plotter_h
#define Plotter_h

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





class Plotter{
    
public:
    
    /// Constructor
    Plotter(const Samples& samples,
            const DrawMode::DrawMode& drawMode);
    
    /// Destructor
    ~Plotter(){};
    
    /// Set options specific for histogram under consideration
    void setOptions(const TString&, const TString&, const TString&, const TString&,
                    const int, const bool, const bool, const bool,
                    const double&, const double&, const double&, const double&,
                    const int, const std::vector<double>&, const std::vector<double>&);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    /// Pair of a legend entry and the histogram for the corresponding sample
    typedef std::pair<TString, TH1*> LegendHistPair;
    
    /// Access histogram under consideration from all samples, do global scalings
    bool prepareDataset(const std::vector<Sample>& v_sample,
                        const std::vector<double>& v_weight);
    
    /// Do stacking, legending, and write in file
    void write(const Channel::Channel& channel, const Systematic::Systematic& systematic);
    
    /// Access global scale factors
    const SystematicChannelFactors& scaleFactors();
    
    
    
    /// Set the style of the plot
    void setStyle(SampleHistPair& sampleHistPair);
    
    /// Set the style of the graph
    void setGraphStyle( TGraph* graph, Style_t marker = 21, Color_t markerColor = 1, Size_t markerSize = 1, 
                        Style_t line = 0, Color_t lineColor = 0, Size_t lineWidth = 1)const;
    
    /// Draw label for decay channel in upper left corner of plot
    void drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize =0.04)const;

    /// Draw official labels (CMS [Preliminary], luminosity and CM energy) above plot
    void drawCmsLabels(const int cmsprelim =1, const double& energy =8, const double& textSize =0.04)const;
    
    /// Draw signal significance label over the plot
    TPaveText* drawSignificance(TH1* signal, TH1* bkg, float Xmin,  float Xmax, float yOffset = 0.f, std::string sLabel ="", const int type=0)const;
    
    /// Calculates purity (type=0) and stability (type=1) curves for a 2D distribution
    TGraphErrors* purityStabilityGraph(TH2* h2d, const int type)const;
    
    /// Calculates binomial uncertainty of the subset/set ratio
    double uncertaintyBinomial(const double pass, const double all)const;
    
    
    
    
    /// Samples to be analysed
    const Samples& samples_;
    
    /// Draw mode for Higgs
    const DrawMode::DrawMode drawMode_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    
    
    /// Pair of histogram to print and corresponding sample
    std::vector<SampleHistPair> v_sampleHistPair_;
    
    
    
    /// Map holding global scale factors for steps they were already accessed, to avoid re-calculating them
    std::map<TString, SystematicChannelFactors> m_stepFactors_;
    
    
    
    /// Name of histogram under consideration
    TString name_;
    TString drawOpt_;
    
    /// Options for the histogram under consideration
    int bins_, rebin_;
    bool stackToNEntries_;
    double rangemin_, rangemax_, ymin_, ymax_;
    std::vector<double> XAxisbins_, XAxisbinCenters_;
    std::vector<Sample::SampleType> sampleTypesToStack2D_;
    TString XAxis_;
    bool logX_, logY_; // The variable logX_ is not used at all...
};







#endif







