#ifndef common_histoUtils_h
#define common_histoUtils_h

#include <vector>

#include <Rtypes.h>

class TH1;
class THStack;
class TStyle;
class TGraphAsymmErrors;
class TGraph;



namespace common{
    
    /// Draw ratio of two histograms
    void drawRatio(const TH1* histNumerator, 
                   const TH1* histDenominator,
                   const TH1 *uncband,
                   const Double_t& ratioMin, 
                   const Double_t& ratioMax, 
                   const bool addFit = 0,
                   const TStyle& myStyle = *gStyle, 
                   const int verbose=0, 
                   const std::vector<double>& err=std::vector<double>(0),
                   const bool useMcStatError = false
                  );

    /// Draw ratio of histograms if not NULL pointers
    void drawRatioXSEC(const TH1* histNumerator, const TH1* histDenominator1, 
                       TGraphAsymmErrors *data_stat = 0, TGraphAsymmErrors *data_syst = 0, 
                       const TH1* histDenominator2 = 0, const TH1* histDenominator3 = 0, const TH1* histDenominator4 = 0, const TH1* histDenominator5 = 0, const TH1* histDenominator6 = 0, const TH1* histDenominator7 = 0,
                       const Double_t& ratioMin = 0.5, const Double_t& ratioMax = 1.5, 
                       TStyle myStyle = *gStyle);//, int verbose = , const std::vector<double>& err=std::vector<double>(0));


    /// Set histogram style to HH definition
    void setHHStyle(TStyle& HHStyle);



    /// Sum the histograms in a stack and return the sum in a new TH1
    TH1* summedStackHisto(const THStack* stack);
    
    /// Add the ratio pad with axis to the specified pad
    TH1* drawRatioPad(TPad* pad, const double yMin, const double yMax, TH1* axisHisto, const double fraction = 0.36, 
                      const TString title = "#frac{Data}{MC}");
    
    /// Get a ratio histogram (errorType: 0-none, 1-nominator, 2-denominator, 3-both)
    TH1* ratioHistogram(const TH1* h_nominator, const TH1* h_denominator, const int errorType = 1);
    
    /// Normalising histogram/graph
    double normalize ( TH1* histo, const double normalization = 1.0, const bool includeOutsideBins = false);
    double normalize ( TGraph* graph, const double normalization = 1.0);
    
    /// Divide each bin by the its width
    void normalizeToBinWidth(TH1* histo);
    
    /// Set the style of the plot
    void setHistoStyle(TH1* hist, Style_t line = 1, Color_t lineColor = 1, Size_t lineWidth = 1, 
                       Style_t fill = 0, Color_t fillColor = 0, 
                       Style_t marker = 0, Color_t markerColor = 1, Size_t markerSize = 1);
    
    /// Set the style of the graph
    void setGraphStyle( TGraph* graph, Style_t marker = 21, Color_t markerColor = 1, Size_t markerSize = 1, 
                        Style_t line = 0, Color_t lineColor = 1, Size_t lineWidth = 1);
}





#endif



