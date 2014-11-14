#include <iostream>
#include <sstream>
#include <cmath>

#include <TH1.h>
#include <TString.h>
#include <THStack.h>
#include <TList.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TF1.h>
#include <TStyle.h>
#include <TText.h>
#include <TGraphAsymmErrors.h>
#include <TExec.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>

#include "plotterUtils.h"





void common::drawRatioXSEC(const TH1* histNumerator, const TH1* histDenominator1, 
                           TGraphAsymmErrors *ratio_stat, TGraphAsymmErrors *ratio_total, 
                           const TH1* histDenominator2, const TH1* histDenominator3, 
                           const TH1* histDenominator4, const TH1* histDenominator5, 
                           const TH1* histDenominator6, const TH1* histDenominator7, 
                           const Double_t& ratioMin, const Double_t& ratioMax, TStyle myStyle)
{
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator_i' (_i = 1-5)
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // the histogram style is transferred from 'histDenominator_i' to the 'ratio_i'
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    /// check that histos exist and have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator1->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        return;
    }

    /// create ratio
    TH1F *ratio1 = 0, *ratio2 = 0, *ratio3 = 0, *ratio4 = 0, *ratio5 = 0, *ratio6 = 0, *ratio7 = 0; 

    ratio1 = (TH1F*)histDenominator1->Clone();
    ratio1->SetLineColor(histDenominator1->GetLineColor());
    ratio1->SetLineStyle(histDenominator1->GetLineStyle());
    ratio1->SetLineWidth(histDenominator1->GetLineWidth());
    ratio1->Divide(histNumerator);

    if (histDenominator2){
        ratio2 = (TH1F*)histDenominator2->Clone();
        ratio2->SetLineColor(histDenominator2->GetLineColor());
        ratio2->SetLineStyle(histDenominator2->GetLineStyle());
        ratio2->SetLineWidth(histDenominator2->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator2->GetNbinsX()){ratio2 = 0;}
        else {ratio2->Divide(histNumerator);}
    };
    if (histDenominator3){
        ratio3 = (TH1F*)histDenominator3->Clone();
        ratio3->SetLineColor(histDenominator3->GetLineColor());
        ratio3->SetLineStyle(histDenominator3->GetLineStyle());
        ratio3->SetLineWidth(histDenominator3->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator3->GetNbinsX()){ratio3 = 0;}
        else {ratio3->Divide(histNumerator);}
    };
    if (histDenominator4){
        ratio4 = (TH1F*)histDenominator4->Clone();
        ratio4->SetLineColor(histDenominator4->GetLineColor());
        ratio4->SetLineStyle(histDenominator4->GetLineStyle());
        ratio4->SetLineWidth(histDenominator4->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator4->GetNbinsX()){ratio4 = 0;}
        else {ratio4->Divide(histNumerator);}
    };
    if (histDenominator5){
        ratio5 = (TH1F*)histDenominator5->Clone();
        ratio5->SetLineColor(histDenominator5->GetLineColor());
        ratio5->SetLineStyle(histDenominator5->GetLineStyle());
        ratio5->SetLineWidth(histDenominator5->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator5->GetNbinsX()){ratio5 = 0;}
        else {ratio5->Divide(histNumerator);}
    };
    if (histDenominator6){
        ratio6 = (TH1F*)histDenominator6->Clone();
        ratio6->SetLineColor(histDenominator6->GetLineColor());
        ratio6->SetLineStyle(histDenominator6->GetLineStyle());
        ratio6->SetLineWidth(histDenominator6->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator6->GetNbinsX()){ratio6 = 0;}
        else {ratio6->Divide(histNumerator);}
    };
    if (histDenominator7){
        ratio7 = (TH1F*)histDenominator7->Clone();
        ratio7->SetLineColor(histDenominator7->GetLineColor());
        ratio7->SetLineStyle(histDenominator7->GetLineStyle());
        ratio7->SetLineWidth(histDenominator7->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator7->GetNbinsX()){ratio7 = 0;}
        else {ratio7->Divide(histNumerator);}
    };

    /// calculate error for ratio only gaussian error of histNumerator
    for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
        if (ratio1) ratio1->SetBinError(bin, sqrt(histDenominator1->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio2) ratio2->SetBinError(bin, sqrt(histDenominator2->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio3) ratio3->SetBinError(bin, sqrt(histDenominator3->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio4) ratio4->SetBinError(bin, sqrt(histDenominator4->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio5) ratio5->SetBinError(bin, sqrt(histDenominator5->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio6) ratio6->SetBinError(bin, sqrt(histDenominator6->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio7) ratio7->SetBinError(bin, sqrt(histDenominator7->GetBinContent(bin))/histNumerator->GetBinContent(bin));
    }

    Int_t    logx  = myStyle.GetOptLogx();
    Double_t left  = myStyle.GetPadLeftMargin();
    Double_t right = myStyle.GetPadRightMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);

    /// create new pad for ratio plot
    TPad *rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    
    /// configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio1->SetStats(kFALSE);
    ratio1->SetTitle("");
    ratio1->SetMaximum(ratioMax);
    ratio1->SetMinimum(ratioMin);
    
    /// configure axis of ratio plot
    ratio1->GetXaxis()->SetTitleSize(histNumerator->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio1->GetXaxis()->SetTitleOffset(histNumerator->GetXaxis()->GetTitleOffset()*0.9);
    ratio1->GetXaxis()->SetLabelSize(histNumerator->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio1->GetXaxis()->SetTitle(histNumerator->GetXaxis()->GetTitle());


    ratio1->GetYaxis()->CenterTitle();
    ratio1->GetYaxis()->SetTitle("#frac{Theory}{Data}");
    ratio1->GetYaxis()->SetTitleSize(histNumerator->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio1->GetYaxis()->SetTitleOffset(histNumerator->GetYaxis()->GetTitleOffset()/scaleFactor);
    ratio1->GetYaxis()->SetLabelSize(histNumerator->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio1->GetYaxis()->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3);
    ratio1->GetYaxis()->SetTickLength(0.03);
    ratio1->GetYaxis()->SetNdivisions(505);
    ratio1->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    ratio1->GetXaxis()->SetNoExponent(kTRUE);
    
    /// delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    histNumerator->GetXaxis()->SetNoExponent(kFALSE);

    /// draw ratio plot
    ratio1->DrawClone("Histo");
    rPad->Update();
    rPad->Modified();
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    
    /// draw grid
    rPad->SetGrid(0,1);

    // draw a horizontal lines on a given histogram
    // a) at 1
    Double_t xmin = ratio1->GetXaxis()->GetXmin();
    Double_t xmax = ratio1->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height.Data(), xmin, xmax);
    f->SetLineStyle(1);//this is frustrating and stupid but apparently necessary...
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    // b) at upper end of ratio pad
    TString height2 = ""; height2 += ratioMax;
    TF1 *f2 = new TF1("f2", height2.Data(), xmin, xmax);
    f2->SetLineStyle(1);
    f2->SetLineWidth(1);
    f2->SetLineColor(kBlack);
    f2->Draw("L same");

    if(ratio_stat) {
        TLegend *leg_band = new TLegend();
        if(ratio_total) leg_band->AddEntry(ratio_total, "Stat. #oplus Syst.", "f");
        leg_band->AddEntry(ratio_stat, "Stat.", "f");
        leg_band->SetX1NDC(0.22);
        leg_band->SetY1NDC(0.97);
        leg_band->SetX2NDC(0.46);
        leg_band->SetY2NDC(0.77);
        leg_band->SetFillStyle(1001);
        leg_band->SetFillColor(10);
        leg_band->SetBorderSize(0);
        leg_band->SetTextSize(0.1);
        leg_band->SetTextAlign(12);
        leg_band->Draw("same");
        if(ratio_total)ratio_total->Draw("same,e2");
        ratio_stat->Draw("same,e2");
    }
    f->Draw("l,same");
    f2->Draw("l,same");
    gPad->RedrawAxis();
    gPad->Update();
    gPad->Modified();
    ratio1->Draw("histo,same");
    if (ratio2) ratio2->Draw("Histo,same");
    if (ratio3) ratio3->Draw("Histo,same");
    if (ratio4) ratio4->Draw("Histo,same");
    if (ratio5) ratio5->Draw("Histo,same");
    if (ratio6) ratio6->Draw("Histo,same");
    if (ratio7) ratio7->Draw("Histo,same");
    gPad->RedrawAxis();
    gPad->Update();
    gPad->Modified();
    rPad->RedrawAxis();
}




void common::drawRatio(const TH1* histNumerator, const TH1* histDenominator, const TH1* uncband,
               const Double_t& ratioMin, const Double_t& ratioMax, 
               bool addFit,
               const TStyle& myStyle, const int verbose, const std::vector<double>& err, const bool useMcStatError)
{
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator'
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = std::sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // if 'err_' is present and its size equals the number of bins in the histos,
    // its valus are considered as error for the ratio
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        return;
    }
    if(verbose>1){
        std::cout << "building ratio plot of " << histNumerator->GetName();
        std::cout << " and " << histDenominator->GetName() << std::endl;
    }

    // create ratio of uncertainty band
    TH1 *band = nullptr;
    if (uncband) {
        band = (TH1*)uncband->Clone("band");
        for(int i=0; i<= 1+uncband->GetNbinsX(); i++)
        {
            double error = 0;
            double content = 1;
            if(band->GetBinContent(i))
            {
                error = band->GetBinError(i) / band->GetBinContent(i);
                content = band->GetBinContent(i) /band->GetBinContent(i);
            }
            band->SetBinError(i, error/content);
            band->SetBinContent(i, content/content);
        }
    }

    // create ratio
    TH1* ratio = (TH1*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    ratio->SetLineColor(histNumerator->GetLineColor());
    ratio->SetLineStyle(histNumerator->GetLineStyle());
    ratio->SetLineWidth(histNumerator->GetLineWidth());
    ratio->SetMarkerColor(histNumerator->GetMarkerColor());
    ratio->SetMarkerStyle(histNumerator->GetMarkerStyle());
    ratio->SetMarkerSize(0.8);
    // calculate error for ratio
    // a) from err_
    if(err.size()==(unsigned int)histNumerator->GetNbinsX()){
        if(verbose>0) std::cout << "ratio error from vector" << std::endl;
        for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
            ratio->SetBinError(bin, err[bin-1]);
        }
    }
    else{
        // b) default: only gaussian error of histNumerator
        if(verbose>0) std::cout << "ratio error from statistical error of " << histNumerator->GetName() << " only" << std::endl;
        for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
            double error = 0.;
            if(!useMcStatError) {
                error = std::sqrt(histNumerator->GetBinContent(bin))/histDenominator->GetBinContent(bin);
            } else {
                double c1 = histNumerator->GetBinContent(bin);
                double c2 = histDenominator->GetBinContent(bin);
                double e1 = histNumerator->GetBinError(bin);
                double e2 = histDenominator->GetBinError(bin);
                error = c1/c2 * sqrt( (e1*e1)/(c1*c1) + (e2*e2)/(c2*c2) );
            }
            ratio->SetBinError(bin, error);
        }
    }
    // get some values from old pad
    //Int_t    logx = gPad->GetLogx();
    //Double_t left = gPad->GetLeftMargin();
    //Double_t right = gPad->GetRightMargin();

    Int_t    logx  = myStyle.GetOptLogx();
    Double_t left  = myStyle.GetPadLeftMargin();
    Double_t right = myStyle.GetPadRightMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);
    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
#ifdef DILEPTON_MACRO
    rPad->SetFillColor(10);
#else
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
#endif
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetName("ratio");
    ratio->SetMaximum(ratioMax);
    ratio->SetMinimum(ratioMin);
    ratio->SetLineWidth(1);
    // configure axis of ratio plot
    ratio->GetXaxis()->SetTitleSize(histNumerator->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio->GetXaxis()->SetTitleOffset(histNumerator->GetXaxis()->GetTitleOffset()*0.9);
    ratio->GetXaxis()->SetLabelSize(histNumerator->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio->GetXaxis()->SetTitle(histNumerator->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetNdivisions(histNumerator->GetNdivisions());
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitle("#frac{N_{Data}}{N_{MC}}");
    ratio->GetYaxis()->SetTitleSize(histNumerator->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio->GetYaxis()->SetTitleOffset(histNumerator->GetYaxis()->GetTitleOffset()/scaleFactor);
    ratio->GetYaxis()->SetLabelSize(histNumerator->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio->GetYaxis()->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3);
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(405);
    ratio->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    // delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    // draw ratio plot
    ratio->DrawClone("p e X0");
    //This is frustrating and stupid but apparently necessary...
    TExec *setex1 { new TExec("setex1","gStyle->SetErrorX(0.5)") };
    setex1->Draw();
    if(band)band->Draw("same,e2");
    TExec *setex2 { new TExec("setex2","gStyle->SetErrorX(0.)") };
    setex2->Draw();
    ratio->SetMarkerSize(0.8);
    ratio->SetLineWidth(2);
    ratio->DrawClone("p e X0 same");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    // draw grid
    rPad->SetGrid(0,1);

    /// make linear fit of the ratio plot
    TF1 *fit = 0;
    if (addFit){
        double ratioxmin = ratio->GetXaxis()->GetXmin();
        double ratioxmax = ratio->GetXaxis()->GetXmax();
        fit = new TF1("fit", "[0]+[1]*x", ratioxmin, ratioxmax);
        fit->SetParName(0, "origin");  fit->SetParameter(0, ratioxmin);
        fit->SetParName(1, "slope");   fit->SetParameter(1, 0);
        fit->SetLineColor(kRed);
        ratio->Fit(fit);
        fit->SetLineColor(kRed);
        fit->Draw("L,same");
        char fitresult[30];
        sprintf(fitresult, "fit = %2.4f + x * %2.4f", fit->GetParameter("origin"), fit->GetParameter("slope"));
        TText *st = new TText (rPad->GetLeftMargin()+0.5, rPad->GetBottomMargin()+0.5, fitresult);
        st->SetTextSize(1.2*st->GetTextSize());
        st->SetNDC(1);
        st->SetTextColor(fit->GetLineColor());
        st->DrawClone();
    }

    // draw a horizontal lines on a given histogram
    // a) at 1
    Double_t xmin = ratio->GetXaxis()->GetXmin();
    Double_t xmax = ratio->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height, xmin, xmax);
    f->SetLineStyle(1);
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    // b) at upper end of ratio pad
    TString height2 = ""; height2 += ratioMax;
    TF1 *f2 = new TF1("f2", height2, xmin, xmax);
    f2->SetLineStyle(1);
    f2->SetLineWidth(1);
    f2->SetLineColor(kBlack);
    f2->Draw("L same");
}



void common::setHHStyle(TStyle& HHStyle, const bool hideErrorX)
{
    const int fontstyle=42;
    HHStyle.SetPalette(1);
        
    // ==============
    //  Canvas
    // ==============
            
    HHStyle.SetCanvasBorderMode(0);
    HHStyle.SetCanvasColor(kWhite);
    HHStyle.SetCanvasDefH(600); //Height of canvas
    HHStyle.SetCanvasDefW(600); //Width of canvas
    HHStyle.SetCanvasDefX(0);   //Position on screen
    HHStyle.SetCanvasDefY(0);
            
    // ==============
    //  Pad
    // ==============
            
    HHStyle.SetPadBorderMode(0);
    // HHStyle.SetPadBorderSize(Width_t size = 1);
    HHStyle.SetPadColor(kWhite);
    HHStyle.SetPadGridX(false);
    HHStyle.SetPadGridY(false);
    HHStyle.SetGridColor(kGray);
    HHStyle.SetGridStyle(3);
    HHStyle.SetGridWidth(1);
            
    // ==============
    //  Frame
    // ==============
            
    HHStyle.SetFrameBorderMode(0);
    HHStyle.SetFrameBorderSize(1);
    HHStyle.SetFrameFillColor(0);
    HHStyle.SetFrameFillStyle(0);
    HHStyle.SetFrameLineColor(1);
    HHStyle.SetFrameLineStyle(1);
    HHStyle.SetFrameLineWidth(1);
            
    // ==============
    //  Histo
    // ==============
    
    if(hideErrorX) HHStyle.SetErrorX(0.0);
    HHStyle.SetEndErrorSize(8);
            
    // HHStyle.SetHistFillColor(1);
    // HHStyle.SetHistFillStyle(0);
    // HHStyle.SetHistLineColor(1);
    HHStyle.SetHistLineStyle(0);
    HHStyle.SetHistLineWidth(1);
    // HHStyle.SetLegoInnerR(Float_t rad = 0.5);
    // HHStyle.SetNumberContours(Int_t number = 20);

    // HHStyle.SetErrorMarker(20);
            
    HHStyle.SetMarkerStyle(20);
            
    // ==============
    //  Fit/function
    // ==============
            
    HHStyle.SetOptFit(1);
    HHStyle.SetFitFormat("5.4g");
    HHStyle.SetFuncColor(2);
    HHStyle.SetFuncStyle(1);
    HHStyle.SetFuncWidth(1);
            
    // ==============
    //  Date
    // ============== 
            
    HHStyle.SetOptDate(0);
    // HHStyle.SetDateX(Float_t x = 0.01);
    // HHStyle.SetDateY(Float_t y = 0.01);
            
    // =====================
    //  Statistics Box
    // =====================
            
    HHStyle.SetOptFile(0);
    HHStyle.SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    HHStyle.SetStatColor(kWhite);
    HHStyle.SetStatFont(fontstyle);
    HHStyle.SetStatFontSize(0.025);
    HHStyle.SetStatTextColor(1);
    HHStyle.SetStatFormat("6.4g");
    HHStyle.SetStatBorderSize(1);
    HHStyle.SetStatH(0.1);
    HHStyle.SetStatW(0.15);
    // HHStyle.SetStatStyle(Style_t style = 1001);
    // HHStyle.SetStatX(Float_t x = 0);
    // HHStyle.SetStatY(Float_t y = 0);
            
    // ==============
    //  Margins
    // ==============

    HHStyle.SetPadTopMargin(0.1);
    HHStyle.SetPadBottomMargin(0.15);
    HHStyle.SetPadLeftMargin(0.20);
    HHStyle.SetPadRightMargin(0.05);
            
    // ==============
    //  Global Title
    // ==============
            
    HHStyle.SetOptTitle(0);
    HHStyle.SetTitleFont(fontstyle);
    HHStyle.SetTitleColor(1);
    HHStyle.SetTitleTextColor(1);
    HHStyle.SetTitleFillColor(10);
    HHStyle.SetTitleFontSize(0.05);
    // HHStyle.SetTitleH(0); // Set the height of the title box
    // HHStyle.SetTitleW(0); // Set the width of the title box
    // HHStyle.SetTitleX(0); // Set the position of the title box
    // HHStyle.SetTitleY(0.985); // Set the position of the title box
    // HHStyle.SetTitleStyle(Style_t style = 1001);
    // HHStyle.SetTitleBorderSize(2);
            
    // ==============
    //  Axis titles
    // ==============
            
    HHStyle.SetTitleColor(1, "XYZ");
    HHStyle.SetTitleFont(fontstyle, "XYZ");
    HHStyle.SetTitleSize(0.05, "XYZ");
    // HHStyle.SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // HHStyle.SetTitleYSize(Float_t size = 0.02);
    HHStyle.SetTitleXOffset(1.0);
    HHStyle.SetTitleYOffset(1.7);
    // HHStyle.SetTitleOffset(1.1, "Y"); // Another way to set the Offset
            
    // ==============
    //  Axis Label
    // ==============
            
    //HHStyle.SetLabelColor(1, "XYZ");
    HHStyle.SetLabelFont(fontstyle, "XYZ");
    HHStyle.SetLabelOffset(0.007, "XYZ");
    HHStyle.SetLabelSize(0.04, "XYZ");
            
    // ==============
    //  Axis
    // ==============
            
    HHStyle.SetAxisColor(1, "XYZ");
    HHStyle.SetStripDecimals(kTRUE);
    HHStyle.SetTickLength(0.03, "XYZ");
    HHStyle.SetNdivisions(510, "XYZ");
    HHStyle.SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    HHStyle.SetPadTickY(1);
            
    // Change for log plots:
    HHStyle.SetOptLogx(0);
    HHStyle.SetOptLogy(0);
    HHStyle.SetOptLogz(0);
            
    // ==============
    //  Text
    // ==============
            
    HHStyle.SetTextAlign(11);
    HHStyle.SetTextAngle(0);
    HHStyle.SetTextColor(1);
    HHStyle.SetTextFont(fontstyle);
    HHStyle.SetTextSize(0.05);
            
    // =====================
    //  Postscript options:
    // =====================
            
    HHStyle.SetPaperSize(20.,20.);
    // HHStyle.SetLineScalePS(Float_t scale = 3);
    // HHStyle.SetLineStyleString(Int_t i, const char* text);
    // HHStyle.SetHeaderPS(const char* header);
    // HHStyle.SetTitlePS(const char* pstitle);
            
    // HHStyle.SetBarOffset(Float_t baroff = 0.5);
    // HHStyle.SetBarWidth(Float_t barwidth = 0.5);
    // HHStyle.SetPaintTextFormat(const char* format = "g");
    // HHStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // HHStyle.SetTimeOffset(Double_t toffset);
    // HHStyle.SetHistMinimumZero(kTRUE);
}



Double_t common::xBinSize(const TAxis* const axis)
{
    // Function taken from http://root.cern.ch/root/html/src/TH1.cxx.html#VwhIWD,
    // but unfortunately not available in ROOT except of in this specific file
    if(!axis->GetNbins()) return 0.;
    bool isEquidistant = true;
    const Double_t firstBinWidth = axis->GetBinWidth(1);
    for(int iBin = 1; iBin < axis->GetNbins()+1; ++iBin){
        const Double_t binWidth = axis->GetBinWidth(iBin);
        const bool match = TMath::AreEqualRel(firstBinWidth, binWidth, TMath::Limits<Double_t>::Epsilon());
        isEquidistant &= match;
        if(!match) break;
    }
    if(isEquidistant) return firstBinWidth;
    else return -999.;
}



void common::addBinDivisionToYaxis(TH1* const hist)
{
    const TAxis* const xAxis = hist->GetXaxis();
    const Double_t binSize = xBinSize(xAxis);
    TString binDivision(" /");
    if(binSize < -990.) binDivision.Append(" #Delta_{bin}");
    else{
        std::ostringstream width;
        width<<" "<<binSize;
        binDivision.Append(width.str());
        TString xTitle = hist->GetXaxis()->GetTitle();
        if(xTitle.EndsWith("]") && xTitle.Contains("[")){
            Ssiz_t last = xTitle.Last('[');
            xTitle = xTitle.Data() + last + 1;
            xTitle.Remove(TString::kTrailing, ']');
            binDivision.Append(" ").Append(xTitle);
        }
        
    }
    TString yTitle = TString(hist->GetYaxis()->GetTitle()).Copy();
    yTitle.Append(binDivision);
    hist->GetYaxis()->SetTitle(yTitle);
}



TH1* common::summedStackHisto(const THStack *stack)
{
    TList* list = stack->GetHists(); //the TList is owned by the stack
    if (list->GetEntries() == 0) return 0;
    TH1* result = (TH1*) list->At(0)->Clone();
    for (int i = 1; i < list->GetEntries(); ++i) {
        result->Add((TH1*)list->At(i));
    }
    return result;
}


TH1* common::drawRatioPad(TPad* pad, const double yMin, const double yMax, TH1* axisHisto, 
                          const TString title, const double fraction)
{
    // y:x size ratio for canvas
    double canvAsym = (pad->GetY2() - pad->GetY1())/(pad->GetX2()-pad->GetX1());
    Double_t left  = pad->GetLeftMargin();
    Double_t right = pad->GetRightMargin();
    // change old pad
    pad->SetBottomMargin(fraction);
    pad->SetRightMargin(right);
    pad->SetLeftMargin(left);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    pad->SetFillColor(10);
    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,fraction+0.001);
#ifdef DILEPTON_MACRO
    rPad->SetFillColor(10);
#else
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
#endif
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*fraction);
    TH1* h_axis = (TH1*)axisHisto->Clone("h_axis");
    h_axis->SetStats(kFALSE);
    h_axis->SetTitle("");
    h_axis->SetName("h_axis");
    h_axis->SetMaximum(yMax);
    h_axis->SetMinimum(yMin);
    // configure axis of the pad
    h_axis->GetXaxis()->SetTitleSize(axisHisto->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    h_axis->GetXaxis()->SetTitleOffset(axisHisto->GetXaxis()->GetTitleOffset()*0.9);
    h_axis->GetXaxis()->SetLabelSize(axisHisto->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    h_axis->GetXaxis()->SetTitle(axisHisto->GetXaxis()->GetTitle());
    h_axis->GetXaxis()->SetNdivisions(axisHisto->GetNdivisions());
    h_axis->GetYaxis()->CenterTitle();
    h_axis->GetYaxis()->SetTitle(title);
    h_axis->GetYaxis()->SetTitleSize(axisHisto->GetYaxis()->GetTitleSize()*scaleFactor);
    h_axis->GetYaxis()->SetTitleOffset(axisHisto->GetYaxis()->GetTitleOffset()/scaleFactor);
    h_axis->GetYaxis()->SetLabelSize(axisHisto->GetYaxis()->GetLabelSize()*scaleFactor);
    h_axis->GetYaxis()->SetLabelOffset(axisHisto->GetYaxis()->GetLabelOffset()*3.3);
    h_axis->GetYaxis()->SetTickLength(0.03);
    h_axis->GetYaxis()->SetNdivisions(405);
    h_axis->GetXaxis()->SetRange(axisHisto->GetXaxis()->GetFirst(), axisHisto->GetXaxis()->GetLast());
    // delete axis of initial plot
    axisHisto->GetXaxis()->SetLabelSize(0);
    axisHisto->GetXaxis()->SetTitleSize(0);
    //This is frustrating and stupid but apparently necessary...
    TExec *setex1 { new TExec("setex1","gStyle->SetErrorX(0.5)") };
    setex1->Draw();
    TExec *setex2 { new TExec("setex2","gStyle->SetErrorX(0.)") };
    setex2->Draw();
    h_axis->Draw("axis");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    pad->SetLeftMargin(left);
    pad->RedrawAxis();
    // draw grid
    rPad->SetGrid(0,1);
    rPad->cd();
    
    // draw a horizontal line on the pad
    Double_t xmin = h_axis->GetXaxis()->GetXmin();
    Double_t xmax = h_axis->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height, xmin, xmax);
    f->SetLineStyle(1);
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    
    
    return h_axis;
    
}


TH1* common::ratioHistogram(const TH1* h_nominator, const TH1* h_denominator, const int errorType)
{
    TH1* h_ratio = (TH1*)h_nominator->Clone();
    for(int iBin = 1; iBin<=h_ratio->GetNbinsX(); ++iBin) {
        const double nominator_v = h_nominator->GetBinContent(iBin);
        const double nominator_e = h_nominator->GetBinError(iBin);
        const double denominator_v = h_denominator->GetBinContent(iBin);
        const double denominator_e = h_denominator->GetBinError(iBin);
        double ratio_e = 0.;
        double ratio_v = denominator_v > 0. ? nominator_v / denominator_v : 0.;
        if(denominator_v > 0.) {
            if(errorType == 1) ratio_e = nominator_e/denominator_v;
            else if(errorType == 2) ratio_e = denominator_e/denominator_v;
            else if(errorType == 3) ratio_e = ratio_v * std::sqrt( (nominator_e*nominator_e)/(nominator_v*nominator_v) + (denominator_e*denominator_e)/(denominator_v*denominator_v) );
        }
        
        h_ratio->SetBinContent(iBin, ratio_v);
        h_ratio->SetBinError(iBin, ratio_e);
    }
    
    return h_ratio;
}


double common::normalize(TH1* histo, const double normalization, const bool includeOutsideBins)
{
    const double integral = includeOutsideBins ? histo->Integral(0, -1) : histo->Integral();
    const double normalizationFactor = normalization / integral;
    histo->Scale(normalizationFactor);
    
    return normalizationFactor;
}


double common::normalize ( TGraph* graph, const double normalization)
{
    const int nBins = graph->GetN();
    
    double integral = 0;
    for(int iBin = 0; iBin < nBins; ++iBin) integral+=graph->GetY()[iBin];
    
    const double normalizationFactor = normalization / integral;
    
    for(int iBin = 0; iBin < nBins; ++iBin) {
        graph->GetY()[iBin] *= normalizationFactor;
        graph->GetEYhigh()[iBin] *= normalizationFactor;
        graph->GetEYlow()[iBin] *= normalizationFactor;
    }
    
    return normalizationFactor;
}


void common::normalizeToBinWidth(TH1* histo)
{
    for(int iBin = 1; iBin <= histo->GetNbinsX(); ++iBin) {
        const double width = histo->GetBinWidth(iBin);
        histo->SetBinContent(iBin, histo->GetBinContent(iBin)/width);
        histo->SetBinError(iBin, histo->GetBinError(iBin)/width);
    }
}


void common::setHistoStyle(TH1* hist, Style_t line, Color_t lineColor, Size_t lineWidth, 
                           Style_t marker, Color_t markerColor, Size_t markerSize,
                           Style_t fill, Color_t fillColor)
{
    if(line != -1) hist->SetLineStyle(line);
    if(lineColor != -1) hist->SetLineColor(lineColor);
    if(lineWidth != -1) hist->SetLineWidth(lineWidth);
    
    if(fill != -1) hist->SetFillStyle(fill);
    if(fillColor != -1) hist->SetFillColor(fillColor);
    
    if(marker != -1) hist->SetMarkerStyle(marker);
    if(markerColor != -1) hist->SetMarkerColor(markerColor);
    if(markerSize != -1) hist->SetMarkerSize(markerSize);
}


void common::setGraphStyle( TGraph* graph, Style_t line, Color_t lineColor, Size_t lineWidth, 
                            Style_t marker, Color_t markerColor, Size_t markerSize,
                            Style_t fill, Color_t fillColor) 
{
    if(line != -1) graph->SetLineStyle(line);
    if(lineColor != -1) graph->SetLineColor(lineColor);
    if(lineWidth != -1) graph->SetLineWidth(lineWidth);
    
    if(fill != -1) graph->SetFillStyle(fill);
    if(fillColor != -1) graph->SetFillColor(fillColor);
    
    if(marker != -1) graph->SetMarkerStyle(marker);
    if(markerColor != -1) graph->SetMarkerColor(markerColor);
    if(markerSize != -1) graph->SetMarkerSize(markerSize);
}


