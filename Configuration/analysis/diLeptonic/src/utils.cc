#include <iostream>
#include <fstream>
#include <sstream>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TPad.h>
#include <TExec.h>
#include <TLine.h>

#include "utils.h"
#include "../../common/include/RootFileReader.h"


void utils::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis)
{
    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}

void utils::drawRatio(TH1* histNumerator, TH1* histDenominator, const TH1* uncband,
               const Double_t& ratioMin, const Double_t& ratioMax,TH1D* hist)
{
    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        std::cout << "building ratio plot of " << histNumerator->GetName();
        std::cout << " and " << histDenominator->GetName() << std::endl;
        return;
    }

    // create ratio of uncertainty band
    TH1 *band = nullptr;
    if (uncband) {
        band = (TH1*)uncband->Clone("band");
        band->Divide(band);
    }


    // create ratio
    TH1* ratio = (TH1*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    ratio->SetLineColor(1);
    ratio->SetLineStyle(histNumerator->GetLineStyle());
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(1);
    // calculate error for ratio
//         for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
//             ratio->SetBinError(bin, std::sqrt(histNumerator->GetBinContent(bin))/histDenominator->GetBinContent(bin));
//         }
    // get some values from old pad
    Int_t    logx = gPad->GetLogx();
    Double_t left = gPad->GetLeftMargin();
    Double_t right = gPad->GetRightMargin();

//     Int_t    logx  = gStyle->GetOptLogx();
//     Double_t left  = gStyle->GetPadLeftMargin();
//     Double_t right = gStyle->GetPadRightMargin();

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
//     // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(1);
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
    //ratio->GetYaxis()->SetTitle("#frac{N_{Data}}{N_{MC}}");
    ratio->GetYaxis()->SetTitle("Data / MC");
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
    ratio->SetMarkerSize(1.0);
    if(hist){
        hist->SetMaximum(ratioMax); 
        hist->SetMinimum(ratioMin); 
        hist->SetTitle(0);
        hist->GetXaxis()->SetLabelSize(0.08);
        hist->GetXaxis()->SetTitleSize(0.08);
        hist->Draw(); 
        for(int i=1;i<hist->GetXaxis()->GetNbins();++i)
        {
            double xLine = hist->GetBinLowEdge(i+1);
            TLine ln(xLine,ratioMin,xLine,ratioMax);
            ln.DrawClone("same");
        }
        ratio->DrawClone("p e X0 same");
    }
    else ratio->DrawClone("p e X0");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    //draw grid
    rPad->SetGrid(0,1);

    // draw a horizontal lines on a given histogram
    // a) at 1
//     Double_t xmin = ratio->GetXaxis()->GetXmin();
//     Double_t xmax = ratio->GetXaxis()->GetXmax();
//     TString height = ""; height += 1;
//     TF1 *f = new TF1("f", height, xmin, xmax);
//     f->SetLineStyle(1);
//     f->SetLineWidth(1);
//     f->SetLineColor(kBlack);
//     f->Draw("L same");
//     // b) at upper end of ratio pad
//     TString height2 = ""; height2 += ratioMax;
//     TF1 *f2 = new TF1("f2", height2, xmin, xmax);
//     f2->SetLineStyle(1);
//     f2->SetLineWidth(1);
//     f2->SetLineColor(kBlack);
//     f2->Draw("L same");
}

void utils::rebin2d(TH2D*& histOut,TH2D* histIn,TString /*name*/,TString xAxisName_,TString yAxisName_,const int rebinX_,const int rebinY_, const int nbinX_, const double* x_binsArr_, const int nbinY_, const double* y_binsArr_)
{
        if(rebinX_||rebinY_){
            histOut = (TH2D* )histIn->Rebin2D(rebinX_,rebinY_,histIn->GetName());
        }
        else
        {
            //histOut = new TH2D(name,"",nbinX_,x_binsArr_,nbinY_,y_binsArr_);
            histOut->GetXaxis()->SetTitle(xAxisName_);
            histOut->GetYaxis()->SetTitle(yAxisName_);
            double binWx = (histIn->GetXaxis())->GetBinWidth(1);
            double binWy = (histIn->GetYaxis())->GetBinWidth(1);
            for(int ix=0;ix<nbinX_;ix++){
                for(int iy=0;iy<nbinY_;iy++){
                    int binx1 = (histIn->GetXaxis())->FindBin(x_binsArr_[ix]+0.1*binWx);
                    int binx2 = (histIn->GetXaxis())->FindBin(x_binsArr_[ix+1]-0.1*binWx);
                    int biny1 = (histIn->GetYaxis())->FindBin(y_binsArr_[iy]+0.1*binWy);
                    int biny2 = (histIn->GetYaxis())->FindBin(y_binsArr_[iy+1]-0.1*binWy);
                    double content=histIn->Integral(binx1,binx2,biny1,biny2);
                    histOut->SetBinContent(ix+1,iy+1,content);
                }
            }
        }
}



TString utils::numToString(double val)
{
    std::stringstream ss;
    ss.str("");
    ss  << val;
    return (TString)ss.str();
}


TString utils::makeBinTitle(TString axisName,double x1,double x2)
{
    std::stringstream ss;
    ss.str("");
    ss  << axisName << " [" << x1 << " : " << x2 << "] ";
    return (TString)ss.str();
    
}

TString utils::makeTitleBins(TString plotNameUnits,std::vector<double>& v_bin,int underflow, int overflow)
{
    std::stringstream ss;
    ss.str("");
    ss  << plotNameUnits << " bins: [";
    if(underflow == 1) ss << "underflow  ";
    for(int i=0; i<(int)v_bin.size(); ++i){
        double bin = v_bin.at(i);
        if(i==0)ss << bin;
        else ss << ",  " << bin;
        //if(i!=0)ss << "  ";
        //ss << "[" << v_bin.at(i) << "  " << v_bin.at(i+1) << "]";
        
    }
    if(overflow == 1) ss << "  overflow";
    ss  << "]";
    return (TString)ss.str();
    
}



void utils::readLineToVector(const TString& file, const TString& keyWord,std::vector<double>& outVector)
{
    std::ifstream tempStream(file.Data(), std::ifstream::in);
    if (!tempStream.good()) {
        std::cerr<<"Error in utils::readLineToVector! Cannot find file with name: "<< file <<"\n...break\n"<<std::endl;
        exit(12);
    }
    while(tempStream.good()){
        std::string line;
        getline(tempStream, line);
        line.erase(0, line.find_first_not_of(" \t"));
        if (line.size() == 0 || line[0] == '#') continue;
        std::vector<TString> vWord;
        std::string word;
        for (std::stringstream ss(line); ss >> word; ){
            vWord.push_back(word);
        }
        if(vWord.at(0) == keyWord.Data())
        {
            vWord.erase(vWord.begin());
            for(auto word: vWord)outVector.push_back(word.Atof());
        }
    }
}



void utils::cat(const TString& file)
{
    std::ifstream tempStream(file.Data(), std::ifstream::in);
    if (!tempStream.good()) {
        std::cerr<<"Error in utils::readLineToVector! Cannot find file with name: "<< file <<"\n...break\n"<<std::endl;
        exit(12);
    }
    while(tempStream.good()){
        std::string line;
        getline(tempStream, line);
        std::cout << line << std::endl;
    }
}



void styleUtils::setResultLegendStyle(TLegend* leg, const bool /*result*/)
{
    double x1 = 0.7, y1 = 0.5;
    double height = 0.2+0.155, width = 0.2;

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void styleUtils::setHHStyle(TStyle& HHStyle)
{
    //const int fontstyle=42;
    HHStyle.SetPalette(1);
        
    // ==============
    //  Canvas
    // ==============
            
//     HHStyle.SetCanvasBorderMode(0);
//     HHStyle.SetCanvasColor(kWhite);
//     HHStyle.SetCanvasDefH(600); //Height of canvas
//     HHStyle.SetCanvasDefW(600); //Width of canvas
//     HHStyle.SetCanvasDefX(0);   //Position on screen
//     HHStyle.SetCanvasDefY(0);
            
    // ==============
    //  Pad
    // ==============
            
//     HHStyle.SetPadBorderMode(0);
//     // HHStyle.SetPadBorderSize(Width_t size = 1);
//     HHStyle.SetPadColor(kWhite);
//     HHStyle.SetPadGridX(false);
//     HHStyle.SetPadGridY(false);
//     HHStyle.SetGridColor(0);
//     HHStyle.SetGridStyle(3);
//     HHStyle.SetGridWidth(1);
            
    // ==============
    //  Frame
    // ==============
            
//     HHStyle.SetFrameBorderMode(0);
//     HHStyle.SetFrameBorderSize(1);
//     HHStyle.SetFrameFillColor(0);
//     HHStyle.SetFrameFillStyle(0);
//     HHStyle.SetFrameLineColor(1);
//     HHStyle.SetFrameLineStyle(1);
//     HHStyle.SetFrameLineWidth(1);
            
    // ==============
    //  Histo
    // ==============

    HHStyle.SetErrorX(0.0);
    HHStyle.SetEndErrorSize(0);
            
    // HHStyle.SetHistFillColor(1);
    // HHStyle.SetHistFillStyle(0);
    HHStyle.SetHistLineColor(1);
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
            
//     HHStyle.SetOptDate(0);
//     // HHStyle.SetDateX(Float_t x = 0.01);
//     // HHStyle.SetDateY(Float_t y = 0.01);
//             
//     // =====================
//     //  Statistics Box
//     // =====================
//             
//     HHStyle.SetOptFile(0);
//     HHStyle.SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
//     HHStyle.SetStatColor(kWhite);
//     HHStyle.SetStatFont(fontstyle);
//     HHStyle.SetStatFontSize(0.025);
//     HHStyle.SetStatTextColor(1);
//     HHStyle.SetStatFormat("6.4g");
//     HHStyle.SetStatBorderSize(1);
//     HHStyle.SetStatH(0.1);
//     HHStyle.SetStatW(0.15);
//     // HHStyle.SetStatStyle(Style_t style = 1001);
//     // HHStyle.SetStatX(Float_t x = 0);
//     // HHStyle.SetStatY(Float_t y = 0);
//             
//     // ==============
//     //  Margins
//     // ==============
// 
//     HHStyle.SetPadTopMargin(0.1);
//     HHStyle.SetPadBottomMargin(0.15);
//     HHStyle.SetPadLeftMargin(0.20);
//     HHStyle.SetPadRightMargin(0.05);
//             
//     // ==============
//     //  Global Title
//     // ==============
//             
//     HHStyle.SetOptTitle(0);
//     HHStyle.SetTitleFont(fontstyle);
//     HHStyle.SetTitleColor(1);
//     HHStyle.SetTitleTextColor(1);
//     HHStyle.SetTitleFillColor(10);
//     HHStyle.SetTitleFontSize(0.05);
//     // HHStyle.SetTitleH(0); // Set the height of the title box
//     // HHStyle.SetTitleW(0); // Set the width of the title box
//     // HHStyle.SetTitleX(0); // Set the position of the title box
//     // HHStyle.SetTitleY(0.985); // Set the position of the title box
//     // HHStyle.SetTitleStyle(Style_t style = 1001);
//     // HHStyle.SetTitleBorderSize(2);
//             
//     // ==============
//     //  Axis titles
//     // ==============
//             
//     HHStyle.SetTitleColor(1, "XYZ");
//     HHStyle.SetTitleFont(fontstyle, "XYZ");
//     HHStyle.SetTitleSize(0.04, "XYZ");
//     // HHStyle.SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
//     // HHStyle.SetTitleYSize(Float_t size = 0.02);
//     HHStyle.SetTitleXOffset(1.25);
//     HHStyle.SetTitleYOffset(1.6);
//     // HHStyle.SetTitleOffset(1.1, "Y"); // Another way to set the Offset
//             
//     // ==============
//     //  Axis Label
//     // ==============
//             
//     //HHStyle.SetLabelColor(1, "XYZ");
//     HHStyle.SetLabelFont(fontstyle, "XYZ");
//     HHStyle.SetLabelOffset(0.007, "XYZ");
//     HHStyle.SetLabelSize(0.04, "XYZ");
//             
//     // ==============
//     //  Axis
//     // ==============
//             
//     HHStyle.SetAxisColor(1, "XYZ");
//     HHStyle.SetStripDecimals(kTRUE);
//     HHStyle.SetTickLength(0.03, "XYZ");
//     HHStyle.SetNdivisions(510, "XYZ");
//     HHStyle.SetPadTickX(1);  // To get tick marks on the opposite side of the frame
//     HHStyle.SetPadTickY(1);
//             
//     // Change for log plots:
//     HHStyle.SetOptLogx(0);
//     HHStyle.SetOptLogy(0);
//     HHStyle.SetOptLogz(0);
//             
//     // ==============
//     //  Text
//     // ==============
//             
//     HHStyle.SetTextAlign(11);
//     HHStyle.SetTextAngle(0);
//     HHStyle.SetTextColor(1);
//     HHStyle.SetTextFont(fontstyle);
//     HHStyle.SetTextSize(0.05);
//             
//     // =====================
//     //  Postscript options:
//     // =====================
//             
//     HHStyle.SetPaperSize(20.,20.);
//     // HHStyle.SetLineScalePS(Float_t scale = 3);
//     // HHStyle.SetLineStyleString(Int_t i, const char* text);
//     // HHStyle.SetHeaderPS(const char* header);
//     // HHStyle.SetTitlePS(const char* pstitle);
//             
//     // HHStyle.SetBarOffset(Float_t baroff = 0.5);
//     // HHStyle.SetBarWidth(Float_t barwidth = 0.5);
//     // HHStyle.SetPaintTextFormat(const char* format = "g");
//     // HHStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
//     // HHStyle.SetTimeOffset(Double_t toffset);
//     // HHStyle.SetHistMinimumZero(kTRUE);
}

