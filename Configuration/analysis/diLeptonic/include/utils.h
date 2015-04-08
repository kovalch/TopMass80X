#ifndef weightUtils_h
#define weightUtils_h

#include <vector>
#include <set>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

namespace styleUtils{
    
    void setHHStyle(TStyle& HHStyle);
    void setResultLegendStyle(TLegend *leg, const bool result);
}


namespace utils{
    
    void addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis);
    void drawRatio(TH1* histNumerator, TH1* histDenominator, const TH1* uncband,const Double_t& ratioMin, const Double_t& ratioMax,TH1D* hist=NULL, TGraphAsymmErrors* graph = NULL ); //TH1D* hist=NULL - needed for all bins in 2d plotting
    void rebin2d(TH2D*& histOut,TH2D* histIn,TString name,TString xAxisName_,TString yAxisName_,const int rebinX_,const int rebinY_, const int nbinX_, const double* x_binsArr_, const int nbinY_, const double* y_binsArr_);
    TString numToString(double val);
    TString makeBinTitle(TString axisName,double x1,double x2);
    TString makeTitleBins(TString plotNameUnits,std::vector<double>& v_bin,int underflow, int overflow);
    
    ///Read line of numbers from txt file and write them to the vector
    void readLineToVector(const TString& file, const TString& keyWord,std::vector<double>& outVector);
    
    ///Printing out whole file
    void cat(const TString& file);
    
    ///Sum of two vectors
    std::vector<double> addVect(const std::vector<double>& a,const std::vector<double>& b, 
                                const double scale=1, const int precision=-999);
    
    ///Difference of two vectors
    std::vector<double> diffVect(const std::vector<double>& a,const std::vector<double>& b, 
                                 const double scale=1, const int precision=-999);
    
    ///Divide vectors
    std::vector<double> divideVect(const std::vector<double>& a,const std::vector<double>& b, 
                                   const double scale=1, const int precision=-999);
    
    ///Relative difference of two vectors 
    std::vector<double> relativeDiff(const std::vector<double>& a,const std::vector<double>& b, 
                                     const double scale=1, const int precision=-999);
    
    ///Set precision for all vector elements
    void setPrecision(std::vector<double>& a,int precision);
    
    /// Prepare canvas and legend
    TCanvas* setCanvas();
    TLegend* setLegend(const double x1 = 0.53,const double y1 =0.60,const double x2 =0.90,const double y2 =0.85);
    
}

#endif

