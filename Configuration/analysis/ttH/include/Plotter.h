#ifndef Plotter_h
#define Plotter_h

#include <vector>
#include <map>

#include <TString.h>

class TLegend;
class RootFileReader;
class TH1;
class TPaveText;

#include "plotterHelpers.h"
#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class AnalysisConfig;





class Plotter{
    
public:
    
    /// Constructor
    Plotter(const char* outputDir,
            const AnalysisConfig& analysisConfig,
            const Samples& samples,
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
    
    /// Draw a label in the top left corner depending on the topLeftLabelId_ parameter
    void drawTopLeftLabel(const TString& text, const double& textSize =0.04)const;
    
    /// Draw signal significance label over the plot
    TPaveText* drawSignificance(const TH1* const signal, const TH1* const signalPlusBackground,
                                const float xMin, const float xMax,
                                const float yOffset = 0.f, const std::string& sLabel ="", const int type =0)const;
                                
    /// Adding histograms to the existing one or creating it
    TH1* addOrCreateHisto(TH1* base, TH1* add_histo);
    
    
    /// Output folder name
    const char* outputDir_;
    
    /// Reference to the analysis config
    const AnalysisConfig& analysisConfig_;
    
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
    
    
    /// Whether pseudodata is used instead of data
    bool hasPseudodata_;
    /// Name of histogram under consideration
    TString name_;
    
    /// Options for the histogram under consideration
    int bins_, rebin_;
    double rangemin_, rangemax_, ymin_, ymax_;
    std::vector<double> XAxisbins_, XAxisbinCenters_;
    TString YAxis_;
    TString XAxis_;
    bool logX_, logY_; // The variable logX_ is not used at all...
    bool normalizeToData_;
};







#endif







