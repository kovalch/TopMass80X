#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>

#include "TString.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/RootFileReader.h"



/// The input base folder
constexpr const char* InputBaseDIR = "selectionRoot";

/// The output base folder
constexpr const char* OutputBaseDIR = "Plots_jetMatch";



void setStyle()
{
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(111110);
    gErrorIgnoreLevel = 1001;
    
    gStyle->SetPalette(1);      //Spektralpalette, Default: 0 resp. 50
    //gStyle->SetNumberContours(20);  // Default: 20
    
    const double width = 600.;
    
    gStyle->SetCanvasDefW(width);
    gStyle->SetCanvasDefH(width);
    
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);
    
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    
    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetTitleOffset(1.2,"X");
    
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    
    // Change for log plots:
    gStyle->SetOptLogx(0);
    gStyle->SetOptLogy(0);
    gStyle->SetOptLogz(0);
    
    TGaxis::SetMaxDigits(3);
    
    //TH1::StatOverflows(kTRUE);// compute mean etc including overflow!
    //gStyle->SetHistMinimumZero(kTRUE); // no zero-suppression on y-axis
    //gStyle->SetOptFit(222);         // 1: Fit-Parameter werden angezeigt
    //gStyle->SetCanvasDefX(400);     // canvas (default) upper edge in X
    //gStyle->SetCanvasDefY(200);     // canvas (default) upper edge in Y
    
    gStyle->SetHistLineWidth(2);
    
    
    //gStyle->SetTitleX(0.2);         // move upper left corner of title box to specified value
    //gStyle->SetTitleY(0.99);        // move upper left corner of title box to specified value
    
    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetLabelSize(0.05,"XY");
}



//Plot one TH1D histogram on each canvas
void plot_TH1(RootFileReader* rootFileReader,
              const TString& outputFolder, const TString& fileName,
              const TString& legendName, const TString& plotName,
              const TString& legendEntry,
              const bool logY =false)
{
    TH1* hist = rootFileReader->GetClone<TH1D>(fileName, plotName);
    
    // Set up histogram style  
    setStyle();
      
    // Set up Canvas
    TCanvas* canvas(0);
    canvas = new TCanvas("canvas");
    canvas->cd();
    
    // Set up Legend
    TLegend* legend(0);
    legend = new TLegend(0.196309, 0.905923, 0.501678, 0.996516);
    legend->SetLineColor(kWhite); 
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);
    legend->SetMargin(0.12);
    legend->SetHeader(legendName);
    legend->AddEntry(hist, legendEntry, "l");
    
    // Scale histograms and get minimum and maximum value 
    //hist1->Scale(1./hist1->Integral(0, hist1->GetNbinsX()+1));
    const Double_t yMax = hist->GetBinContent(hist->GetMaximumBin());
    if(logY){
        hist->GetYaxis()->SetRangeUser(0.1, 4.*yMax);
        canvas->SetLogy();
    }
    else{
        hist->GetYaxis()->SetRangeUser(0., 1.4*yMax);
    }
    
    // Draw the histograms
    hist->SetLineWidth(2);
    hist->SetTitle("");
    hist->Draw();
            
    canvas->Modified();
    canvas->Update();

    // Adjust stats boxes, has to be done only after drawing
    TPaveStats* stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    const double x2ndc = 0.99;
    stats->SetX2NDC(x2ndc);
    stats->SetX1NDC(x2ndc-0.2);
    //const double y2ndc = 0.99;
    //stats->SetY2NDC(y2ndc);
    //stats->SetY1NDC(y2ndc-0.2);
    stats->SetLineColor(hist->GetLineColor());
    
    canvas->Modified();
    canvas->Update();
    
    legend->Draw("same"); 
    canvas->Modified();
    canvas->Update();
    
    canvas->Print(outputFolder+plotName+".eps");
    
    legend->Delete();
    canvas->Close();
}

   
   
//Plot one TH2D histogram on each canvas   
void plot_TH2(RootFileReader* rootFileReader,
              const TString& outputFolder, const TString& fileName,
              const TString& legendName, const TString& plotName,
              const TString& legendEntry)
{
    TH1* hist = rootFileReader->GetClone<TH2D>(fileName, plotName);
    
    // Set up histogram style  
    setStyle();
    gStyle->SetOptStat("e");
    
    // Set up canvas
    TCanvas* canvas(0);
    canvas = new TCanvas("canvas");
    canvas->cd();
    
    hist->Draw("COLZ");
    
    TLegend* legend(0);
    legend = new TLegend(0.196309, 0.905923, 0.501678, 0.996516);
    legend->SetLineColor(kWhite); 
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);
    legend->SetHeader(legendName);
    legend->AddEntry(hist, legendEntry, "l");
    legend->Draw("same"); 
    
    canvas->Print(outputFolder+plotName+".eps");
    
    legend->Delete();
    canvas->Close();
}


    
void printHistogram(const TString& fileName, const TString& legendName, const TString& step,
                    const std::vector<Channel::Channel>& v_channel,
                    const std::vector<Systematic::Systematic>& v_systematic)
{
    TString name;
    
    const TString prefix = "jetMatch_";
    
    const std::vector<TString> v_whichSelection = {
        "initial_",
        "initialAmbiguous_",
        "initialUnambiguous_",
        "mismatchedInR_",
        "matchedInR_",
        "matchedInRAmbiguous_",
        "matchedInRUnambiguous_",
        
        "mismatchedInPt_m05_p07_",
        "mismatchedInPtAmbiguousMatchedInR_m05_p07_",
        "mismatchedInPtUnambiguousMatchedInR_m05_p07_",
        "mismatchedInPtAmbiguous_m05_p07_",
        "mismatchedInPtUnambiguous_m05_p07_",
        "matchedInPt_m05_p07_",
        "matchedInPtAmbiguous_m05_p07_",
        "matchedInPtUnambiguous_m05_p07_",
        
        "mismatchedInPt_m05_p06_",
        "mismatchedInPtAmbiguousMatchedInR_m05_p06_",
        "mismatchedInPtUnambiguousMatchedInR_m05_p06_",
        "mismatchedInPtAmbiguous_m05_p06_",
        "mismatchedInPtUnambiguous_m05_p06_",
        "matchedInPt_m05_p06_",
        "matchedInPtAmbiguous_m05_p06_",
        "matchedInPtUnambiguous_m05_p06_",
        
        "mismatchedInPt_m05_p05_",
        "mismatchedInPtAmbiguousMatchedInR_m05_p05_",
        "mismatchedInPtUnambiguousMatchedInR_m05_p05_",
        "mismatchedInPtAmbiguous_m05_p05_",
        "mismatchedInPtUnambiguous_m05_p05_",
        "matchedInPt_m05_p05_",
        "matchedInPtAmbiguous_m05_p05_",
        "matchedInPtUnambiguous_m05_p05_",
        
        "mismatchedInPt_m04_p07_",
        "mismatchedInPtAmbiguousMatchedInR_m04_p07_",
        "mismatchedInPtUnambiguousMatchedInR_m04_p07_",
        "mismatchedInPtAmbiguous_m04_p07_",
        "mismatchedInPtUnambiguous_m04_p07_",
        "matchedInPt_m04_p07_",
        "matchedInPtAmbiguous_m04_p07_",
        "matchedInPtUnambiguous_m04_p07_",
        
        "mismatchedInPt_m04_p06_",
        "mismatchedInPtAmbiguousMatchedInR_m04_p06_",
        "mismatchedInPtUnambiguousMatchedInR_m04_p06_",
        "mismatchedInPtAmbiguous_m04_p06_",
        "mismatchedInPtUnambiguous_m04_p06_",
        "matchedInPt_m04_p06_",
        "matchedInPtAmbiguous_m04_p06_",
        "matchedInPtUnambiguous_m04_p06_",
        
        "mismatchedInPt_m04_p05_",
        "mismatchedInPtAmbiguousMatchedInR_m04_p05_",
        "mismatchedInPtUnambiguousMatchedInR_m04_p05_",
        "mismatchedInPtAmbiguous_m04_p05_",
        "mismatchedInPtUnambiguous_m04_p05_",
        "matchedInPt_m04_p05_",
        "matchedInPtAmbiguous_m04_p05_",
        "matchedInPtUnambiguous_m04_p05_",
        
        "mismatchedInPt_m03_p07_",
        "mismatchedInPtAmbiguousMatchedInR_m03_p07_",
        "mismatchedInPtUnambiguousMatchedInR_m03_p07_",
        "mismatchedInPtAmbiguous_m03_p07_",
        "mismatchedInPtUnambiguous_m03_p07_",
        "matchedInPt_m03_p07_",
        "matchedInPtAmbiguous_m03_p07_",
        "matchedInPtUnambiguous_m03_p07_",
        
        "mismatchedInPt_m03_p06_",
        "mismatchedInPtAmbiguousMatchedInR_m03_p06_",
        "mismatchedInPtUnambiguousMatchedInR_m03_p06_",
        "mismatchedInPtAmbiguous_m03_p06_",
        "mismatchedInPtUnambiguous_m03_p06_",
        "matchedInPt_m03_p06_",
        "matchedInPtAmbiguous_m03_p06_",
        "matchedInPtUnambiguous_m03_p06_",
        
        "mismatchedInPt_m03_p05_",
        "mismatchedInPtAmbiguousMatchedInR_m03_p05_",
        "mismatchedInPtUnambiguousMatchedInR_m03_p05_",
        "mismatchedInPtAmbiguous_m03_p05_",
        "mismatchedInPtUnambiguous_m03_p05_",
        "matchedInPt_m03_p05_",
        "matchedInPtAmbiguous_m03_p05_",
        "matchedInPtUnambiguous_m03_p05_",
        
        
    };
    
    const std::vector<TString> v_whichJets = {
        "all_",
        "topORhiggs_",
        "top_",
        "higgs_",
    };
    
    RootFileReader* rootFileReader = RootFileReader::getInstance();
    
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
            const TString inputFolder = common::accessFolder(InputBaseDIR, channel, systematic);
            const TString outputFolder = common::assignFolder(OutputBaseDIR, channel, systematic);
            TString path = outputFolder;
            path.Append(fileName).Append("/");
            gSystem->MakeDirectory(path.Data());
            
            const TString channelName = Channel::convert(channel);
            const TString inputFilename = inputFolder + channelName + "_" + fileName + ".root";
            
            for(const TString& whichSelection : v_whichSelection){
                for(const TString& whichJets : v_whichJets){
                    const TString identifier = whichSelection + whichJets;
                    
                    name = identifier + "deltaR";
                    plot_TH1(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets, true);
                    
                    name = identifier + "deltaR_zoom";
                    plot_TH1(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets, true);
                    
                    name = identifier + "recoPtOverGenPt";
                    plot_TH1(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets, true);
                    
                    name = identifier + "deltaPtOverGenPt";
                    plot_TH1(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets, true);
                    
                    name = identifier + "genPtVsRecoPt";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "genEtaVsRecoEta";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "genPhiVsRecoPhi";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "genPVsRecoP";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "genPtVsDeltaPt";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "deltaRVsGenPt";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "recoPtOverGenPtVsGenPt";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "recoPtOverGenPtVsGenEta";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "recoPtOverGenPtVsGenPhi";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "recoPtOverGenPtVsDeltaR";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                    name = identifier + "deltaPtOverGenPtVsDeltaR";
                    plot_TH2(rootFileReader, path, inputFilename, channelName+" "+legendName, prefix+name+step, whichJets);
                    
                }
            }
        }
    }
}



/// All systematics allowed as steering parameter for executable
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal,
    };
}

int main(int argc, char** argv){
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convert(Channel::allowedChannelsPlotting)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal", false, 1, 100,
        common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLAnalyser::interpretGlobal(argc, argv);
        
    //Hardcoded input files, since detailed validation makes sense only on samples containing ttbar system, and some only on real ttH events
    const std::vector<TString> v_inputFileTtbar = {
        //"ttbarsignalPlusBbbar",
        //"ttbarsignalPlusOther",
        "ttbarH125tobbbar"
    };
    
    const std::vector<TString> v_inputLegTtbar = {
        //"t#bar{t}signal+b#bar{b}",
        //"t#bar{t}signal+Other",
        "t#bar{t}H_{125}->b#bar{b}"
    };
    
    const std::vector<TString> v_step = {
        "_step7"
    };  
    
    // Set up channels
    //std::vector<Channel::Channel> v_channel(Channel::allowedChannelsPlotting);
    std::vector<Channel::Channel> v_channel = {Channel::mumu, Channel::emu, Channel::ee};
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for(auto i_channel : v_channel) std::cout << Channel::convert(i_channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
    if(opt_systematic.isSet()) v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    std::cout << "Processing systematics: ";
    for (auto systematic : v_systematic) std::cout << systematic.name() << " ";
    std::cout << "\n\n";
    
    // Loop over all samples and plot the histograms
    for(size_t iFileName = 0; iFileName < v_inputFileTtbar.size(); ++iFileName){
        for(const auto& step : v_step){
            const TString& fileName = v_inputFileTtbar.at(iFileName);
            const TString& legendName = v_inputLegTtbar.at(iFileName);

            printHistogram(fileName, legendName, step, v_channel, v_systematic);
        }
    }
}
