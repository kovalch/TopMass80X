#ifndef common_histoUtils_h
#define common_histoUtils_h

#include <vector>

#include <Rtypes.h>

class TH1;
class THStack;
class TStyle;
class TGraphAsymmErrors;



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
                   const std::vector<double>& err=std::vector<double>(0));

    /// Draw ratio of histograms if not NULL pointers
    void drawRatioXSEC(const TH1* histNumerator, const TH1* histDenominator1, 
                       TGraphAsymmErrors *data_stat = 0, TGraphAsymmErrors *data_syst = 0, 
                       const TH1* histDenominator2 = 0, const TH1* histDenominator3 = 0, const TH1* histDenominator4 = 0, const TH1* histDenominator5 = 0, const TH1* histDenominator6 = 0, const TH1* histDenominator7 = 0,
                       const Double_t& ratioMin = 0.5, const Double_t& ratioMax = 1.5, 
                       TStyle myStyle = *gStyle);//, int verbose = , const std::vector<double>& err=std::vector<double>(0));


    /// Set histogram style to HH definition
    void setHHStyle(TStyle& HHStyle);



    /** Sum the histograms in a stack and
     * return the sum in a new TH1
     */
    TH1* summedStackHisto(const THStack* stack);
}





#endif



