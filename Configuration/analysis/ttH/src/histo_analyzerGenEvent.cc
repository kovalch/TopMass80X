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
constexpr const char* OutputBaseDIR = "Plots_genEvent";



void setStyle()
{
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptTitle(0);
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
    legend = new TLegend(0.5759, 0.669774, 0.882119, 0.799263);
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



//Plot two TH1D histograms on the same canvas               
void plot_2TH1(RootFileReader* rootFileReader,
               const TString& outputFolder, const TString& fileName,
               const TString& legendName, const TString& plotName1, const TString& plotName2,
               const TString& legendEntry1, const TString& legendEntry2,
               const bool logY =false)
{
    TH1* hist1 = rootFileReader->GetClone<TH1D>(fileName, plotName1);
    TH1* hist2 = rootFileReader->GetClone<TH1D>(fileName, plotName2);
    
    // Set up histogram style  
    setStyle();
    
    // Set up Canvas
    TCanvas* canvas(0);
    canvas = new TCanvas("canvas");
    canvas->cd();
    
    // Set up Legend
    TLegend* legend(0);
    legend = new TLegend(0.5759, 0.669774, 0.882119, 0.799263);
    legend->SetLineColor(kWhite); 
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);
    legend->SetMargin(0.12);
    legend->SetHeader(legendName);
    legend->AddEntry(hist1, legendEntry1, "l");
    legend->AddEntry(hist2, legendEntry2, "l");
    
    // set up histograms
    hist1->SetLineColor(1);
    hist2->SetLineColor(2);
    std::vector<TH1*> v_hist;
    v_hist.push_back(hist1);
    v_hist.push_back(hist2);
    
    // Scale histograms and get minimum and maximum value
    Double_t yMin(99999.);
    Double_t yMax(-99999.);
    for(size_t i = 0; i != v_hist.size(); ++i){    
        //v_hist.at(i)->Scale(1./v_hist.at(i)->Integral(0, v_hist.at(i)->GetNbinsX()+1));
        const Double_t yMinHist = v_hist.at(i)->GetBinContent(v_hist.at(i)->GetMinimumBin());
        const Double_t yMaxHist = v_hist.at(i)->GetBinContent(v_hist.at(i)->GetMaximumBin());
        if(yMin > yMinHist) yMin = yMinHist;
        if(yMax < yMaxHist) yMax = yMaxHist;
    }
    
    
    // Draw the histograms
    for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
        TH1* hist = *i_hist;
        const bool firstHist = i_hist==v_hist.begin();
        hist->SetLineWidth(2);
        
        if(firstHist){
            if(logY){
                hist->GetYaxis()->SetRangeUser(0.1, 4.*yMax);
                canvas->SetLogy();
            }
            else{
                hist->GetYaxis()->SetRangeUser(0., 1.4*yMax);
            }
            
            hist->SetTitle("");
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
    
    TString plotName = plotName1;
    plotName.ReplaceAll("1st", "");
    canvas->Print(outputFolder+plotName+".eps");

    legend->Delete();
    canvas->Close();
}



//Plot four TH1D histograms on the same canvas
void plot_4TH1(RootFileReader* rootFileReader,
               const TString& outputFolder, const TString& fileName,
               const TString& legendName, const TString& plotName1, const TString& plotName2, const TString& plotName3, const TString& plotName4,
               const TString& legendEntry1, const TString& legendEntry2, const TString& legendEntry3, const TString& legendEntry4,
               const bool logY =false)
{
    TH1* hist1 = rootFileReader->GetClone<TH1D>(fileName, plotName1);
    TH1* hist2 = rootFileReader->GetClone<TH1D>(fileName, plotName2);
    TH1* hist3 = rootFileReader->GetClone<TH1D>(fileName, plotName3);
    TH1* hist4 = rootFileReader->GetClone<TH1D>(fileName, plotName4);
    
    // Set up histogram style  
    setStyle();
    
    // Set up Canvas
    TCanvas* canvas(0);
    canvas = new TCanvas("canvas");
    canvas->cd();
    
    // Set up Legend
    TLegend* legend(0);
    legend = new TLegend(0.610802, 0.590036, 0.914863, 0.775565);
    legend->SetLineColor(kWhite);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.038);
    legend->SetMargin(0.12);
    legend->SetHeader(legendName);
    legend->AddEntry(hist1, legendEntry1, "l");
    legend->AddEntry(hist2, legendEntry2, "l");
    legend->AddEntry(hist3, legendEntry3, "l");
    legend->AddEntry(hist4, legendEntry4, "l");
    
    // set up histograms
    hist1->SetLineColor(1);
    hist2->SetLineColor(2);
    hist3->SetLineColor(3);
    hist4->SetLineColor(4);
    std::vector<TH1*> v_hist;
    v_hist.push_back(hist1);
    v_hist.push_back(hist2);
    v_hist.push_back(hist3);
    v_hist.push_back(hist4);
    
    // Scale histograms and get minimum and maximum value
    Double_t yMin(99999.);
    Double_t yMax(-99999.);
    for(size_t i=0 ; i!=v_hist.size(); i++){    
        //v_hist.at(i)->Scale(1./v_hist.at(i)->Integral(0, v_hist.at(i)->GetNbinsX()+1));
        const Double_t yMinHist = v_hist.at(i)->GetBinContent(v_hist.at(i)->GetMinimumBin());
        const Double_t yMaxHist = v_hist.at(i)->GetBinContent(v_hist.at(i)->GetMaximumBin());
        if(yMin > yMinHist) yMin = yMinHist;
        if(yMax < yMaxHist) yMax = yMaxHist;
    }
    
    // Draw the histograms
    for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
        TH1* hist = *i_hist;
        const bool firstHist = i_hist==v_hist.begin();
        hist->SetLineWidth(2);
        
        if(firstHist){
            if(logY){
                hist->GetYaxis()->SetRangeUser(0.1, 4.*yMax);
                canvas->SetLogy();
            }
            else{
                hist->GetYaxis()->SetRangeUser(0., 1.4*yMax);
            }
            
            hist->SetTitle("");
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
    
    TString plotName = plotName1;
    plotName.ReplaceAll("1st", "");
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
    legend = new TLegend(0.5759, 0.669774, 0.882119, 0.799263);
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
    TString name, name1, name2, name3, name4;
    
    const TString prefix = "genEvent_";
    
    const std::vector<TString> v_topOrHiggsName = {
        "top_",
        "higgs_",
    };
    
    const std::vector<TString> v_topOrHiggsNameForBjets = {
        "allB_",
        "top_",
        "higgs_",
        "topHiggs_",
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
            
            name = "top_matchedBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "top_");
            
            name = "top_unmatchedGenBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "top_");
            
            name = "top_unmatchedRecoBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "top_");
            
            name = "higgs_matchedBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "higgs_");
            
            name = "higgs_unmatchedGenBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "higgs_");
            
            name = "higgs_unmatchedRecoBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "higgs_");
            
            name = "topHiggs_matchedBjet";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "topHiggs_");
            
            
            for(const TString& topOrHiggsName : v_topOrHiggsNameForBjets){
                 name = topOrHiggsName + "bhadrons_multiplicity";
                 plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                 
                 name = topOrHiggsName + "bjets_multiplicity";
                 plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                 
                 name = topOrHiggsName + "bhadronsPerBjet_multiplicity";
                 plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                 
                 name = topOrHiggsName + "bhadronsVsBjets_multiplicity";
                 plot_TH2(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
            }
            
            for(const TString& topOrHiggsName : v_topOrHiggsName){
                name = topOrHiggsName + "dijet_deltaEta";
                plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                
                name = topOrHiggsName + "dijet_deltaPhi";
                plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                
                name = topOrHiggsName + "dijet_deltaR";
                plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, topOrHiggsName);
                
                name1 = topOrHiggsName + "jet1st_pt";
                name2 = topOrHiggsName + "jet2nd_pt";
                plot_2TH1(rootFileReader, path, inputFilename, legendName, prefix+name1+step, prefix+name2+step, topOrHiggsName+"leading", topOrHiggsName+"subleading");
                
                name1 = topOrHiggsName + "jet1st_eta";
                name2 = topOrHiggsName + "jet2nd_eta";
                plot_2TH1(rootFileReader, path, inputFilename, legendName, prefix+name1+step, prefix+name2+step, topOrHiggsName+"leading", topOrHiggsName+"subleading");
            }
            
            name = "topHiggs_dijet_deltaRMin";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "topHiggs_");
            
            name = "topHiggs_dijet_deltaEtaMin";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "topHiggs_");
            
            name = "topHiggs_dijet_deltaPhiMin";
            plot_TH1(rootFileReader, path, inputFilename, legendName, prefix+name+step, "topHiggs_");
            
            name1 = "topHiggs_topJet1st_pt";
            name2 = "topHiggs_topJet2nd_pt";
            name3 = "topHiggs_higgsJet1st_pt";
            name4 = "topHiggs_higgsJet2nd_pt";
            plot_4TH1(rootFileReader, path, inputFilename, legendName, prefix+name1+step,  prefix+name2+step, prefix+name3+step, prefix+name4+step, "top_leading", "top_subleading", "higgs_leading", "higgs_subleading");
            
            name1 = "topHiggs_topJet1st_eta";
            name2 = "topHiggs_topJet2nd_eta";
            name3 = "topHiggs_higgsJet1st_eta";
            name4 = "topHiggs_higgsJet2nd_eta";
            plot_4TH1(rootFileReader, path, inputFilename, legendName, prefix+name1+step,  prefix+name2+step, prefix+name3+step, prefix+name4+step, "top_leading", "top_subleading", "higgs_leading", "higgs_subleading");
        }
    }
}



int main(int argc, char** argv){
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convert(Channel::allowedChannelsPlotting)));    
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic({Systematic::nominal});
    
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
    std::vector<Channel::Channel> v_channel = {Channel::mumu, Channel::emu, Channel::ee};
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for(auto i_channel : v_channel) std::cout << Channel::convert(i_channel) << " ";
    std::cout << "\n\n";
    
    // Loop over all samples and plot the histograms
    for(size_t iFileName = 0; iFileName < v_inputFileTtbar.size(); ++iFileName){
        for(const auto& step : v_step){
            const TString& fileName = v_inputFileTtbar.at(iFileName);
            const TString& legendName = v_inputLegTtbar.at(iFileName);

            printHistogram(fileName, legendName , step, v_channel, v_systematic);
        }
    } 
}




