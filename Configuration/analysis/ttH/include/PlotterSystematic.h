#ifndef PlotterSystematic_h
#define PlotterSystematic_h

#include <vector>
#include <set>
#include <map>

#include <TString.h>

class TLegend;
class RootFileReader;
class TH1;
class TH2;
class TGraph;
class TGraphErrors;
class TPaveText;
class TPad;

#include "plotterHelpers.h"
#include "SamplesFwd.h"
#include "Sample.h"
#include "../../common/include/sampleHelpers.h"





class PlotterSystematic{
    
public:
    
    /// Constructor
    PlotterSystematic(const char* outputDir, 
                      const std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > >& inputFileLists);
    
    /// Destructor
    ~PlotterSystematic(){};
    
    /// Set options specific for histogram under consideration
    void setOptions(const TString&, const TString&, const TString&, const TString&,
                    const int, const bool, const bool, const bool,
                    const double&, const double&, const double&, const double&);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    /// Pair of a legend entry and the histogram for the corresponding sample
    typedef std::pair<TH1*, TH1*> HistoPair;
    typedef std::map<Systematic::Type, HistoPair> SystematicHistoMap;
    
    /// Prepare styling parameters
    void prepareStyle();

    /// Calculate the differential cross section
    void writeVariations(const SystematicHistoMap& histoCollection, const Channel::Channel channel, const TString processName, const bool logY = false);
    
    /// Set the style of the plot
    void setHistoStyle(TH1* hist, Style_t line = 1, Color_t lineColor = 1, Size_t lineWidth = 1, 
                       Style_t fill = 0, Color_t fillColor = 0, 
                       Style_t marker = 0, Color_t markerColor = 1, Size_t markerSize = 1)const;
    
    /// Set the style of the graph
    void setGraphStyle( TGraph* graph, Style_t marker = 21, Color_t markerColor = 1, Size_t markerSize = 1, 
                        Style_t line = 0, Color_t lineColor = 1, Size_t lineWidth = 1)const;

    /// Update histogram axis
    void updateHistoAxis(TH1* histo, const bool logY=false)const;
    
    /// Draw label for decay channel in upper left corner of plot
    void drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize =0.04)const;

    /// Draw official labels (CMS [Preliminary], luminosity and CM energy) above plot
    void drawCmsLabels(const int cmsprelim =1, const double& energy =8, const double& textSize =0.04)const;
    
    /// Normalise histogram to the unit area
    void normalize(TH1* histo)const;
    
    /// Add area for the ratio plot
    TH1* drawRatioPad(TPad* pad, const double yMin, const double yMax, TH1* axisHisto, const double fraction = 0.36, 
                      const TString title = "#frac{Systematic}{Nominal}")const;
    
    /// Get a ratio histogram
    TH1* ratioHistogram(const TH1* h_nominator, const TH1* h_denominator, const int errorId = 0)const;
    
    /// Calculates binomial uncertainty of the subset/set ratio
    double uncertaintyBinomial(const double pass, const double all)const;
    
    

    /// Output folder name
    const char* outputDir_;
    
    /// Samples to be analysed
    const std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > >& inputFileLists_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    /// Styles of lines for UP/DOWN systematic variations
    std::pair<Style_t, Style_t> lineStyles_;
    
    /// Colors for systematic variations
    std::vector<Color_t> lineColors_;
    
    /// Name of histogram under consideration
    TString name_;
    
    /// Name of process under consideration
    TString processName_;
    
    /// Vector of process names that should be plotted
    std::vector<TString> processNames_;
    
    
    /// Options for the histogram under consideration
    int signalType_;
    bool plotResponse_;
    double rangemin_, rangemax_, ymin_, ymax_;
//     std::vector<double> XAxisbins_, XAxisbinCenters_;

    TString YAxis_, XAxis_;
    bool logX_, logY_; // The variable logX_ is not used at all...
    
    /// Data for cross section calculation
    struct ValueError {
        double v; double e;
        ValueError(double val, double err):v(val), e(err){}
        ValueError():v(1.), e(1.){}
        double v2()const {return v*v;}
        double e2()const {return e*e;}
        double eOv()const {return e/v;}
        double eOv2()const {return (e*e)/(v*v);}
    };
};







#endif







