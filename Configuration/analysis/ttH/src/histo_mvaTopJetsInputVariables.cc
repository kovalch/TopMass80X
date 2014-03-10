#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"

#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/utils.h"
#include "../../common/include/RootFileReader.h"



/// The input base folder
constexpr const char* InputBaseDIR = "mvaOutput";

/// The output base folder
constexpr const char* OutputBaseDIR = "Plots_mvaVariables";





void setStyle(){
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptStat(111110);
    gErrorIgnoreLevel = 1001;

    gStyle->SetPalette(1);      //Spektralpalette, Default: 0 resp. 50
    //gStyle->SetNumberContours(20);  // Default: 20

    double width = 600.;

    gStyle->SetCanvasDefW(width);
    gStyle->SetCanvasDefH(width);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetTitleOffset(1.2,"X");


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



void printHistogram(const TString& inputFile,
                    const TString& outputFile,
                    const std::pair<TString, TString>& nameStepPair,
                    const std::vector<std::string>& v_drawMode)
{
    const TString& histnameCorrect = nameStepPair.first;
    TString histnameSwapped(histnameCorrect);
    histnameSwapped.ReplaceAll("correct_", "swapped_");
    TString histnameWrong(histnameCorrect);
    histnameWrong.ReplaceAll("correct_", "wrong_");
    
    // Set up canvas
    TCanvas* canvas(0);
    canvas = new TCanvas("canvas1");
    canvas->cd();
    
    // Set up legend
    TLegend* legend(0);
    legend = new TLegend(0.16, 0.84, 0.36, 0.99);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->SetMargin(0.12);
    legend->SetFillStyle(1001);
    
    // set up histograms
    std::vector<TH1*> v_hist;
    {
        RootFileReader* rootFileReader(RootFileReader::getInstance());
        TH1* histCorrect = rootFileReader->GetClone<TH1D>(inputFile, histnameCorrect);
        histCorrect->SetLineColor(4);
        legend->AddEntry(histCorrect, "correct", "l");
        v_hist.push_back(histCorrect);
        if(std::find(v_drawMode.begin(), v_drawMode.end(), "swap") != v_drawMode.end()){
            TH1* histSwapped = rootFileReader->GetClone<TH1D>(inputFile, histnameSwapped);
            histSwapped->SetLineColor(3);
            legend->AddEntry(histSwapped, "swapped", "l");
            v_hist.push_back(histSwapped);
        }
        TH1* histWrong = rootFileReader->GetClone<TH1D>(inputFile, histnameWrong);
        histWrong->SetLineColor(2);
        legend->AddEntry(histWrong, "wrong", "l");
        v_hist.push_back(histWrong);
    }
    
    // Scale histograms and get minimum and maximum value
    Double_t yMin(99999.);
    Double_t yMax(-99999.);
    for(TH1* hist : v_hist){
        hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
        const Double_t yMinHist = hist->GetMinimum();
        const Double_t yMaxHist = hist->GetMaximum();
        if(yMin > yMinHist) yMin = yMinHist;
        if(yMax < yMaxHist) yMax = yMaxHist;
    }
    
    // Draw the histograms
    for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
        TH1* hist = *i_hist;
        const bool firstHist = i_hist==v_hist.begin();
        hist->SetLineWidth(2);
        
        if(firstHist){
            hist->SetTitle("");
            hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
            hist->Draw();
            
            canvas->Modified();
            canvas->Update();
        }
        else{
            hist->Draw("sameS");
            
            canvas->Modified();
            canvas->Update();
        }
    }
    
    
    // Adjust stats boxes, has to be done only after drawing
    std::vector<TPaveStats*> v_stats;
    for(size_t iHist = 0; iHist < v_hist.size(); ++iHist){
        TH1* hist = v_hist.at(iHist);
        
        TPaveStats* stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
        const double x2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        stats->SetX2NDC(x2ndc);
        stats->SetX1NDC(x2ndc-0.2);
        //const double y2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        //stats->SetY2NDC(y2ndc);
        //stats->SetY1NDC(y2ndc-0.2);
        stats->SetLineColor(hist->GetLineColor());
        
        canvas->Modified();
        canvas->Update();
    }
    
    legend->Draw("same");
    canvas->Modified();
    canvas->Update();
    
    const TString plotName(outputFile);
    
    canvas->Print(plotName + ".eps");
    canvas->Print(plotName + ".png");
    
    legend->Delete();
    for(size_t iHist = 0; iHist < v_hist.size(); ++iHist) v_hist.at(iHist)->Delete();
    v_hist.clear();
    for(size_t iStats = 0; iStats < v_stats.size(); ++iStats) v_stats.at(iStats)->Delete();
    v_stats.clear();
    canvas->Close();
    
    
}



void histoMvaTopJetsinputVariables(const std::vector<Channel::Channel>& v_channel,
                                   const std::vector<Systematic::Systematic>& v_systematic,
                                   const std::vector<std::string>& v_drawMode)
{
    
    constexpr const char* prefix = "mvaP_";
    
    // Histograms to be printed
    const std::vector<TString> v_histonameBase = {
        "jetChargeDiff",
        "meanDeltaPhi_b_met",
        "massDiff_recoil_bbbar",
        "pt_b_antiLepton",
        "pt_antiB_lepton",
        "deltaR_b_antiLepton",
        "deltaR_antiB_lepton",
        "btagDiscriminatorSum",
        "deltaPhi_antiBLepton_bAntiLepton",
        "massDiff_fullBLepton_bbbar",
        "meanMt_b_met",
        "massSum_antiBLepton_bAntiLepton",
        "massDiff_antiBLepton_bAntiLepton",
    };
    
    // Set up histogram style
    setStyle();
    
    // Loop over all systematics and channels
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
            
            // Compose input filename
            const TString inputFolder = common::accessFolder(InputBaseDIR, channel, systematic);
            const TString inputFile = inputFolder + "plots.root";
            
            // Compose output folder
            const TString outputFolder = common::assignFolder(OutputBaseDIR, channel, systematic);
            
            // Loop over all histogams
            for(const TString& histonameBase : v_histonameBase){
                
                // Find histogram names and selection steps from histogram for correct combinations
                TString histoname("correct_");
                histoname.Append(prefix).Append(histonameBase);
                std::cout<<"Producing histograms for: "<<histonameBase<<"\n";
                const std::vector<std::pair<TString, TString> > v_nameStepPair =
                    tth::nameStepPairs(inputFile, histoname);
                std::cout<<std::endl;
                
                // Loop over steps and print histograms
                for(const auto& nameStepPair : v_nameStepPair){
                    const TString outputFile = outputFolder + histonameBase + nameStepPair.second;
                    printHistogram(inputFile, outputFile, nameStepPair, v_drawMode);
                }
            }
        }
    }
    
    std::cout<<"MVA input variables plotting successfully finished\n";
}



int main(int argc, char** argv)
{
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convertChannels(Channel::allowedChannelsPlotting)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal", false, 1, 100,
        common::makeStringCheck({"Nominal", ""}));
    CLParameter<std::string> opt_drawMode("d", "Draw modes - whether to draw curve for swapped combinations (swap)", false, 1, 1,
        common::makeStringCheck({"swap", ""}));
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up channels
    std::vector<Channel::Channel> v_channel({Channel::combined});
    if(opt_channel.isSet()) v_channel = Channel::convertChannels(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for (auto channel : v_channel)std::cout << Channel::convertChannel(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic({Systematic::nominal});
    if(opt_systematic.isSet()) v_systematic = Systematic::convertSystematics(opt_systematic.getArguments());
    std::cout << "Processing systematics: ";
    for (auto systematic : v_systematic) std::cout << Systematic::convertSystematic(systematic) << " ";
    std::cout << "\n\n";
    
    // Set up draw modes
    std::vector<std::string> v_drawMode({""});
    if(opt_drawMode.isSet()) v_drawMode = opt_drawMode.getArguments();
    std::cout << "Processing draw modes: ";
    for (auto mode : v_drawMode) std::cout << mode << " ";
    std::cout << "\n\n";
    
    histoMvaTopJetsinputVariables(v_channel, v_systematic, v_drawMode);
}





