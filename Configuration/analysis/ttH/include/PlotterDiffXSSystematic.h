#ifndef PlotterDiffXSSystematic_h
#define PlotterDiffXSSystematic_h

#include <vector>
#include <set>
#include <map>

#include <TString.h>

class RootFileReader;
class TH1;
class TH2;
class TGraphAsymmErrors;
class TPad;

#include "plotterHelpers.h"
#include "SamplesFwd.h"
#include "Sample.h"
#include "../../common/include/sampleHelpers.h"





class PlotterDiffXSSystematic{
    
public:
    
    /// Constructor
    PlotterDiffXSSystematic(const char* outputDir, 
                      const std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > >& inputFileLists);
    
    /// Destructor
    ~PlotterDiffXSSystematic(){};
    
    /// Set options specific for histogram under consideration
    void setOptions(const TString&, const TString&, const TString&, const TString&,
                    const int, const bool, const bool, const bool,
                    const double&, const double&, const double&, const double&);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    enum ErrorType{stat, syst, total, shape, rate};
        
    /// A pair of up/down variations (relative)
    struct UpDown {
        double u; double d;
        UpDown(double up, double down):u(up), d(down){}
        UpDown():u(0.), d(0.){}
        
        void square() {
            u = u*u;
            d = d*d;
        }
        
        void sqrt() {
            u = std::sqrt(u);
            d = std::sqrt(d);
        }
        
        void addInQuadrature(UpDown err) {
            u += err.u*err.u;
            d += err.d*err.d;
        }
        
        double maxAbsVariation()const {
            return std::max(std::fabs(u), std::fabs(d));
        }
    };
    
    // A list of systematic types that should be combined into a single systematic
    // Combination types: 0 - added in quadrature; 1 - largest absolute variation in each direction;
    struct SystematicCombination {
        int type; std::vector<Systematic::Type> systematics;
        SystematicCombination():type(1.), systematics(0){};
        SystematicCombination(int combinationType):type(combinationType), systematics(0){};
        SystematicCombination(int combinationType, Systematic::Type systematicTypes[]):type(combinationType){
            systematics = std::vector<Systematic::Type>( systematicTypes, systematicTypes + sizeof(systematicTypes) / sizeof(systematicTypes[0]) );
        };
        
        void addSystematic(Systematic::Type systematicType) {
            systematics.push_back(systematicType);
        }
        
        bool hasSystematic(Systematic::Type systematicType) {
            return std::find(systematics.begin(), systematics.end(), systematicType) != systematics.end();
        }
    };
    
    /// List of up/down errors for each uncertainty type
    typedef std::map<Systematic::Type, UpDown> ErrorMap;
    
    /// Pair of a legend entry and the histogram for the corresponding sample
    typedef std::pair<TH1*, TH1*> HistoPair;
    
    /// List of histograms for each systematic
    typedef std::map<Systematic::Type, HistoPair> SystematicHistoMap;
    
//     /// Pair of systematic type and systematics-combination type
//     typedef std::pair<Systematic::Type, int> SystematicCombination;
    
    /// Prepare styling parameters
    void prepareStyle();
    
    /// Read hisotgrams with a particular name from all systematic source files
    SystematicHistoMap readSystematicHistograms(TString histoName, const Channel::Channel& channel)const;
    
    /// Get single up/down (1/-1) histogram out of many different PDF variations
    TH1* getPdfHisto(TString histoName, const Channel::Channel& channel, const int variation)const;

    /// Plot different systematic shapes for each process
    std::vector<ErrorMap> extractVariations(const SystematicHistoMap& m_systematicHistos)const;
    
    /// Plot the final cross section with uncertainties
    void plotXSection(const Channel::Channel& channel);
    
    /// Set errors to each of the histogram
    TGraphAsymmErrors* errorGraphFromHisto(const TH1* histo, const std::vector<UpDown>& errors)const;

    /// Update histogram axis
    void updateHistoAxis(TPad* pad)const;
    
    /// Calculates uncertainty of particular type for a single bin
    UpDown binUncertaintyOfType(const ErrorMap& m_errors, const ErrorType type = ErrorType::total)const;
    
    /// Calculated uncertainties from each type of variation for a single bin
    ErrorMap binUncertaintiesOfType(const ErrorMap& m_errors, const ErrorType type = ErrorType::total)const;
    
    /// Print individual uncertainties for each bin and mean uncertainty across all bins
    void printAllUncertainties(const std::vector<ErrorMap>& errorMaps, const ErrorMap& errorMap_inclusive)const;
    
    /// Calculates binomial uncertainty of the subset/set ratio
    double uncertaintyBinomial(const double pass, const double all)const;
    
    /// Function for sorting uncertainties that are used to calculate a median: largest absolute values first
    static bool uncertaintySortingFunction(const double a, const double b);
    
    

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
    
    /// Vector of process names that should be plotted
    std::vector<Systematic::Type> normalisedSystematicTypes_;
    
    /// A pair of [new systematic, combination type] for each systematic that sould be merged with others
    std::map<Systematic::Type, SystematicCombination> combinedSystematicTypes_;
    
    
    /// Options for the histogram under consideration
    double rangemin_, rangemax_, ymin_, ymax_;
//     std::vector<double> XAxisbins_, XAxisbinCenters_;
    int nRatio_max_;

    TString YAxis_, XAxis_;
    bool logX_, logY_; // The variable logX_ is not used at all...
};







#endif







