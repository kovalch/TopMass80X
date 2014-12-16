#ifndef weightUtils_h
#define weightUtils_h

#include <vector>
#include <set>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>


namespace styleUtils{
    
    void setHHStyle(TStyle& HHStyle);
    void setResultLegendStyle(TLegend *leg, const bool result);
}


namespace utils{
    
    void addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis);
    void drawRatio(TH1* histNumerator, TH1* histDenominator, const TH1* uncband,const Double_t& ratioMin, const Double_t& ratioMax,TH1D* hist=NULL);
    void rebin2d(TH2D*& histOut,TH2D* histIn,TString name,TString xAxisName_,TString yAxisName_,const int rebinX_,const int rebinY_, const int nbinX_, const double* x_binsArr_, const int nbinY_, const double* y_binsArr_);
    TString numToString(double val);
    TString makeBinTitle(TString axisName,double x1,double x2);
    TString makeTitleBins(TString plotNameUnits,std::vector<double>& v_bin,int underflow, int overflow);
}




#endif

