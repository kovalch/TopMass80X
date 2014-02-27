#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "TString.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"

/// The input base folder
constexpr const char* InputBaseDIR = "selectionRoot/Nominal/";

/// The output base folder
constexpr const char* OutputBaseDIR = "Plots_genLevel_genRecoMatching";


void setStyle(){
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
void plot_TH1 (TString outputFolder, TH1D *hist1, TString fileName, TString LegName, TString FigName, 
               int hist1_color, TString hist1_leg, bool LogY, bool zoomXaxis=false, double xl=0.0, double xu=0.0) {
    
    // Set up histogram style  
    setStyle();
      
    // Set up Canvas
    TCanvas* C1(0);
    C1 = new TCanvas("C1","C1",400,400);
    C1->Divide(1,1);
    C1->cd(1);
    
    if(LogY) gPad->SetLogy();
    
    // Set up Legend
    TLegend* leg(0);
    leg = new TLegend(0.5759,0.669774,0.882119,0.799263);
    leg->SetLineColor(kWhite); 
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetMargin(0.12);
    leg->SetHeader(LegName);
    
    // set up histograms
    hist1->SetLineColor(hist1_color);
    leg->AddEntry(hist1, hist1_leg, "l");
    
    // Scale histograms and get minimum and maximum value 
    //hist1->Scale(1./hist1->Integral(0, hist1->GetNbinsX()+1));
    const Double_t yMax = hist1->GetBinContent(hist1->GetMaximumBin());
    
       
    // Draw the histograms
    hist1->SetLineWidth(2);
    hist1->SetTitle("");
    if (!LogY) hist1->GetYaxis()->SetRangeUser(0.,1.4*yMax);
    if (zoomXaxis) hist1->GetXaxis()->SetRangeUser(xl,xu);
    hist1->Draw();
            
    C1->Modified();
    C1->Update();


    // Adjust stats boxes, has to be done only after drawing
    std::vector<TPaveStats*> v_stats;

    TPaveStats* stats = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");
        const double x2ndc = 0.99;
        stats->SetX2NDC(x2ndc);
        stats->SetX1NDC(x2ndc-0.2);
        //const double y2ndc = 0.99;
        //stats->SetY2NDC(y2ndc);
        //stats->SetY1NDC(y2ndc-0.2);
        stats->SetLineColor(hist1->GetLineColor());
        
        C1->Modified();
        C1->Update();
    
    
    leg->Draw("same"); 
    C1->Modified();
    C1->Update();
    
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".eps");
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".png");

    leg->Delete();
    C1->Close();
                     
}

//Plot two TH1D histograms on the same canvas               
void plot_2TH1 (TString outputFolder, TH1D *hist1, TH1D *hist2, TString fileName, TString LegName, TString FigName, 
                int hist1_color, int hist2_color, TString hist1_leg, TString hist2_leg, 
                bool zoomXaxis=false, double xl=0.0, double xu=0.0) {

    // Set up histogram style  
    setStyle();
        
    // Set up Canvas
    TCanvas* C1(0);
    C1 = new TCanvas("C1","C1",400,400);
    C1->Divide(1,1);
    C1->cd(1);
            
    // Set up Legend
    TLegend* leg(0);
    leg = new TLegend(0.5759,0.669774,0.882119,0.799263);
    leg->SetLineColor(kWhite); 
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetMargin(0.12);
    leg->SetHeader(LegName);
            
    // set up histograms
    hist1->SetLineColor(hist1_color);
    hist2->SetLineColor(hist2_color);
    leg->AddEntry(hist1, hist1_leg, "l");
    leg->AddEntry(hist2, hist2_leg, "l");
    
    std::vector<TH1D*> v_hist;{
            v_hist.push_back(hist1);
            v_hist.push_back(hist2);
    }
    
    // Scale histograms and get minimum and maximum value
    Double_t yMin(99999.);
    Double_t yMax(-99999.);
    for(size_t i=0 ; i!=v_hist.size(); i++){    
        //v_hist.at(i)->Scale(1./v_hist.at(i)->Integral(0, v_hist.at(i)->GetNbinsX()+1));
        const Double_t yMinHist = v_hist.at(i)->GetMinimum();
        const Double_t yMaxHist = v_hist.at(i)->GetMaximum();
        if(yMin > yMinHist) yMin = yMinHist;
        if(yMax < yMaxHist) yMax = yMaxHist;
     }
     
          
    // Draw the histograms
    for(std::vector<TH1D*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
        TH1D* hist = *i_hist;
        const bool firstHist = i_hist==v_hist.begin();
        hist->SetLineWidth(2);
        
        if(firstHist){
            hist->SetTitle("");
            hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
            if (zoomXaxis) hist->GetXaxis()->SetRangeUser(xl,xu);
            hist->Draw();
            
            C1->Modified();
            C1->Update();
        }
        else{
            hist->Draw("sameS");
                    
            C1->Modified();
            C1->Update();
        }
    }
     
    // Adjust stats boxes, has to be done only after drawing
    std::vector<TPaveStats*> v_stats;
    for(size_t iHist = 0; iHist < v_hist.size(); ++iHist){
        TH1D* hist = v_hist.at(iHist);
                 
        TPaveStats* stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
        const double x2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        stats->SetX2NDC(x2ndc);
        stats->SetX1NDC(x2ndc-0.2);
        //const double y2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        //stats->SetY2NDC(y2ndc);
        //stats->SetY1NDC(y2ndc-0.2);
        stats->SetLineColor(hist->GetLineColor());
                 
        C1->Modified();
        C1->Update();
    }
             
    leg->Draw("same"); 
    C1->Modified();
    C1->Update();
    
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".eps");
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".png");

    leg->Delete();
    C1->Close();

}


//Plot four TH1D histograms on the same canvas
void plot_4TH1 (TString outputFolder, TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4, TString fileName, TString LegName, TString FigName, 
                int hist1_color, int hist2_color, int hist3_color, int hist4_color,
                TString hist1_leg, TString hist2_leg, TString hist3_leg, TString hist4_leg,
                bool zoomXaxis=false, double xl=0.0, double xu=0.0) {    
    // Set up histogram style  
    setStyle();
   
    // Set up Canvas
    TCanvas* C1(0);
    C1 = new TCanvas("C1","C1",400,400);
    C1->Divide(1,1);
    C1->cd(1);
   
    // Set up Legend
    TLegend* leg(0);
    leg = new TLegend(0.610802,0.590036,0.914863,0.775565);
    leg->SetLineColor(kWhite);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.038);
    leg->SetMargin(0.12);
    leg->SetHeader(LegName);
 
    // set up histograms
    hist1->SetLineColor(hist1_color);
    hist2->SetLineColor(hist2_color);
    hist3->SetLineColor(hist3_color);
    hist4->SetLineColor(hist4_color);
    leg->AddEntry(hist1, hist1_leg, "l");
    leg->AddEntry(hist2, hist2_leg, "l");
    leg->AddEntry(hist3, hist3_leg, "l");
    leg->AddEntry(hist4, hist4_leg, "l");
    std::vector<TH1D*> v_hist;{
        v_hist.push_back(hist1);
        v_hist.push_back(hist2);
        v_hist.push_back(hist3);
        v_hist.push_back(hist4);
    }

    // Scale histograms and get minimum and maximum value
    Double_t yMin(99999.);
    Double_t yMax(-99999.);
    for(size_t i=0 ; i!=v_hist.size(); i++){    
    //v_hist.at(i)->Scale(1./v_hist.at(i)->Integral(0, v_hist.at(i)->GetNbinsX()+1));
    const Double_t yMinHist = v_hist.at(i)->GetMinimum();
    const Double_t yMaxHist = v_hist.at(i)->GetMaximum();
    if(yMin > yMinHist) yMin = yMinHist;
    if(yMax < yMaxHist) yMax = yMaxHist;
    }

 
    // Draw the histograms
    for(std::vector<TH1D*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
        TH1D* hist = *i_hist;
        const bool firstHist = i_hist==v_hist.begin();
        hist->SetLineWidth(2);

        if(firstHist){
            hist->SetTitle("");
            hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
            if (zoomXaxis) hist->GetXaxis()->SetRangeUser(xl,xu);
            hist->Draw();
   
            C1->Modified();
            C1->Update();
        }
        else{
            hist->Draw("sameS");

            C1->Modified();
            C1->Update();
        }
    }
   
    // Adjust stats boxes, has to be done only after drawing
    std::vector<TPaveStats*> v_stats;
    for(size_t iHist = 0; iHist < v_hist.size(); ++iHist){
        TH1D* hist = v_hist.at(iHist);
   
        TPaveStats* stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
        const double x2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        stats->SetX2NDC(x2ndc);
        stats->SetX1NDC(x2ndc-0.2);
        //const double y2ndc = 0.99 - static_cast<double>(iHist)*0.21;
        //stats->SetY2NDC(y2ndc);
        //stats->SetY1NDC(y2ndc-0.2);
        stats->SetLineColor(hist->GetLineColor());
     
        C1->Modified();
        C1->Update();
    }

    leg->Draw("same"); 
    C1->Modified();
    C1->Update();
    
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".eps");
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".png");

    leg->Delete();
    C1->Close();
    
}
   
   
//Plot one TH2D histogram on each canvas   
void plot_TH2 (TString outputFolder, TH2D *hist1, TString fileName, TString LegName, TString FigName, TString hist1_leg, 
               bool zoomXaxis=false, double xl=0.0, double xu=0.0, bool zoomYaxis=false, double yl=0.0, double yu=0.0) {
    // Set up histogram style  
    setStyle();
    gStyle->SetOptStat("e");
    TCanvas* C1(0);
      
    C1 = new TCanvas("C1","C1",400,400);
    C1->Divide(1,1);
    C1->cd(1);
    if (zoomXaxis) hist1->GetXaxis()->SetRangeUser(xl,xu);
    if (zoomYaxis) hist1->GetYaxis()->SetRangeUser(yl,yu);
    hist1->Draw("COLZ");
    
    TLegend* leg(0);
    leg = new TLegend(0.5759,0.669774,0.882119,0.799263);
    leg->SetLineColor(kWhite); 
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetHeader(LegName);
    leg->AddEntry(hist1, hist1_leg, "l");
    leg->Draw("same"); 
    
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".eps");
    C1->SaveAs(outputFolder+fileName+"_"+FigName+".png");

    leg->Delete();
    C1->Close();
    
}
    
void printHistogram(TString fileName, TString legFileName, TString prefix_histname, TString ending_histname, bool mismatched, const std::vector<Channel::Channel>& v_channel, const std::vector<Systematic::Systematic>& v_systematic){
    
    for(const auto& i_systematic : v_systematic){
        for(const auto& i_channel : v_channel){
            const TString outputFolder = common::assignFolder(OutputBaseDIR, i_channel, i_systematic);
            TString channel          = Channel::convertChannel(i_channel);
            
            TString name=channel+"_"+fileName+ending_histname;
            TString name_leg=channel+" "+legFileName;
            TString file_name= InputBaseDIR+channel+"/"+channel+"_"+fileName+".root";
            
            
            //Open the file with the histograms
            TFile *file = TFile::Open(file_name);
            
            //If you want to print the histograms for the Gen-Reco Matching where at least 2 Gen Jets are matched to the same Reco Jets
            //then you have to change the "name" as following:
            if (mismatched) name=channel+"_"+fileName+"_Mismatched";
        
            if (!mismatched){
                //set the names for the histograms
                std::vector<TString> v_Overlapping_All_histname1D_;{
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"matchedBjetFromTop"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"unmatchedGenBjetFromTop"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"unmatchedRecoBjetFromTop"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"matchedBjetFromHiggs"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"unmatchedGenBjetFromHiggs"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"unmatchedRecoBjetFromHiggs"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"matchedBjet"+ending_histname);
                }

                // set up histograms
                std::vector<TH1D*> v_Overlapping_All_histname1D;
                
                for(size_t iHist = 0; iHist < v_Overlapping_All_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Overlapping_All_histname1D_.at(iHist));
                    v_Overlapping_All_histname1D.push_back(histname1D);
                }
                
                //Print Histograms
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(0), name, name_leg, "matchedBjetFromTop", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "unmatchedGenBjetFromTop", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "unmatchedRecoBjetFromTop", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "matchedBjetFromHiggs", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "unmatchedGenBjetFromHiggs", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "unmatchedRecoBjetFromHiggs", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "matchedBjet", 8, "Top/Higgs Jets", false);
                
            }
            if (!mismatched){
                
                //=============================        Plots for testing overlapping   ================================//
                
                //-----------Plots for all Jets (No matter if they are Top or Higgs)-----------//
                //set the names for the histograms
                std::vector<TString> v_Overlapping_All_histname1D_;{
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"HasOverlappingHadrons"+ending_histname);
                    v_Overlapping_All_histname1D_.push_back(prefix_histname+"NOverlappingHadrons"+ending_histname);
                }
                std::vector<TString> v_Overlapping_All_histname2D_;{
                    v_Overlapping_All_histname2D_.push_back(prefix_histname+"GenJetsID_VS_NOverlappingHadrons_2D"+ending_histname);
                }
                
                // set up histograms
                std::vector<TH1D*> v_Overlapping_All_histname1D;
                std::vector<TH2D*> v_Overlapping_All_histname2D;
                
                for(size_t iHist = 0; iHist < v_Overlapping_All_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Overlapping_All_histname1D_.at(iHist));
                    v_Overlapping_All_histname1D.push_back(histname1D);
                }
                for(size_t iHist = 0; iHist < v_Overlapping_All_histname2D_.size(); ++iHist){
                    TH2D *histname2D  = (TH2D*) file->Get(v_Overlapping_All_histname2D_.at(iHist));
                    v_Overlapping_All_histname2D.push_back(histname2D);
                }
 
                //Print Histograms
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(0), name, name_leg, "HasOverlappingHadrons_All", 8, "All Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_All_histname1D.at(1), name, name_leg, "NOverlappingHadrons_All", 8, "All Jets", false);
                plot_TH2 (outputFolder, v_Overlapping_All_histname2D.at(0), name, name_leg, "GenJetsID_VS_NOverlappingHadrons_2D_All","All Jets");
           
                //----------------            Plots for the Top Jets           -------------------//
                //set the names for the histograms
                std::vector<TString> v_Overlapping_Top_histname1D_;{
                    v_Overlapping_Top_histname1D_.push_back(prefix_histname+"HasOverlappingHadrons_Top"+ending_histname);
                    v_Overlapping_Top_histname1D_.push_back(prefix_histname+"NOverlappingHadrons_Top"+ending_histname);
                }
                std::vector<TString> v_Overlapping_Top_histname2D_;{
                    v_Overlapping_Top_histname2D_.push_back(prefix_histname+"GenJetsID_VS_NOverlappingHadrons_Top_2D"+ending_histname);
                }

                // set up histograms
                std::vector<TH1D*> v_Overlapping_Top_histname1D;
                std::vector<TH2D*> v_Overlapping_Top_histname2D;
                for(size_t iHist = 0; iHist < v_Overlapping_Top_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Overlapping_Top_histname1D_.at(iHist));
                    v_Overlapping_Top_histname1D.push_back(histname1D);
                }
                for(size_t iHist = 0; iHist < v_Overlapping_Top_histname2D_.size(); ++iHist){
                    TH2D *histname2D  = (TH2D*) file->Get(v_Overlapping_Top_histname2D_.at(iHist));
                    v_Overlapping_Top_histname2D.push_back(histname2D);
                } 
                
                //Print Histograms
                plot_TH1 (outputFolder, v_Overlapping_Top_histname1D.at(0), name, name_leg, "HasOverlappingHadrons_Top", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_Top_histname1D.at(1), name, name_leg, "NOverlappingHadrons_Top", 4, "Top Jets", false);
                plot_TH2 (outputFolder, v_Overlapping_Top_histname2D.at(0), name, name_leg, "GenJetsID_VS_NOverlappingHadrons_2D_Top","Top Jets");
  
                
                //-----------------           Plots for the Higgs Jets         ---------------------//
                //set the names for the histograms
                std::vector<TString> v_Overlapping_Higgs_histname1D_;{
                    v_Overlapping_Higgs_histname1D_.push_back(prefix_histname+"HasOverlappingHadrons_Higgs"+ending_histname);
                    v_Overlapping_Higgs_histname1D_.push_back(prefix_histname+"NOverlappingHadrons_Higgs"+ending_histname);
                }
                std::vector<TString> v_Overlapping_Higgs_histname2D_;{
                    v_Overlapping_Higgs_histname2D_.push_back(prefix_histname+"GenJetsID_VS_NOverlappingHadrons_Higgs_2D"+ending_histname);
                }
   
                // set up histograms
                std::vector<TH1D*> v_Overlapping_Higgs_histname1D;
                std::vector<TH2D*> v_Overlapping_Higgs_histname2D;
                for(size_t iHist = 0; iHist < v_Overlapping_Higgs_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Overlapping_Higgs_histname1D_.at(iHist));
                    v_Overlapping_Higgs_histname1D.push_back(histname1D);
                }
                for(size_t iHist = 0; iHist < v_Overlapping_Higgs_histname2D_.size(); ++iHist){
                    TH2D *histname2D  = (TH2D*) file->Get(v_Overlapping_Higgs_histname2D_.at(iHist));
                    v_Overlapping_Higgs_histname2D.push_back(histname2D);
                }  
                
                //Print Histograms
                plot_TH1 (outputFolder, v_Overlapping_Higgs_histname1D.at(0), name, name_leg, "HasOverlappingHadrons_Higgs", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_Higgs_histname1D.at(1), name, name_leg, "NOverlappingHadrons_Higgs", 7, "Higgs Jets", false);
                plot_TH2 (outputFolder, v_Overlapping_Higgs_histname2D.at(0), name, name_leg, "GenJetsID_VS_NOverlappingHadrons_2D_Higgs","Higgs Jets");

                
                //-------------------        Plots for Top + Higgs Jets   ---------------------//
                //set the names for the histograms
                std::vector<TString> v_Overlapping_TopHiggs_histname1D_;{
                    v_Overlapping_TopHiggs_histname1D_.push_back(prefix_histname+"HasOverlappingHadrons_TopHiggs"+ending_histname);
                    v_Overlapping_TopHiggs_histname1D_.push_back(prefix_histname+"NOverlappingHadrons_TopHiggs"+ending_histname);
                }
                std::vector<TString> v_Overlapping_TopHiggs_histname2D_;{
                    v_Overlapping_TopHiggs_histname2D_.push_back(prefix_histname+"GenJetsID_VS_NOverlappingHadrons_TopHiggs_2D"+ending_histname);
                }
                
                // set up histograms
                std::vector<TH1D*> v_Overlapping_TopHiggs_histname1D;
                std::vector<TH2D*> v_Overlapping_TopHiggs_histname2D;
                for(size_t iHist = 0; iHist < v_Overlapping_TopHiggs_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Overlapping_TopHiggs_histname1D_.at(iHist));
                    v_Overlapping_TopHiggs_histname1D.push_back(histname1D);
                }
                for(size_t iHist = 0; iHist < v_Overlapping_TopHiggs_histname2D_.size(); ++iHist){
                    TH2D *histname2D  = (TH2D*) file->Get(v_Overlapping_TopHiggs_histname2D_.at(iHist));
                    v_Overlapping_TopHiggs_histname2D.push_back(histname2D);
                }
                
                //Print Histograms
                plot_TH1 (outputFolder, v_Overlapping_TopHiggs_histname1D.at(0), name, name_leg, "HasOverlappingHadrons_TopHiggs", 8, "Top+Higgs Jets", false);
                plot_TH1 (outputFolder, v_Overlapping_TopHiggs_histname1D.at(1), name, name_leg, "NOverlappingHadrons_TopHiggs", 8, "Top+Higgs Jets", false);
                plot_TH2 (outputFolder, v_Overlapping_TopHiggs_histname2D.at(0), name, name_leg, "GenJetsID_VS_NOverlappingHadrons_2D_TopHiggs","Top+Higgs Jets");
   
                
                //===============                   Gen Level Plots                   =================//
                //----------------------          Plots for the Top Jets        ---------------------------------//
                //set the names for the histograms
                std::vector<TString> v_Gen_Top_histname1D_;{
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"leading_Topjet_Pt"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"subleading_Topjet_Pt"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"leading_Topjet_eta"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"subleading_Topjet_eta"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"DeltaEta_Topjet"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"DeltaPhi_Topjet"+ending_histname);
                    v_Gen_Top_histname1D_.push_back(prefix_histname+"DeltaR_Topjet"+ending_histname);
                }
                
                // set up histograms
                std::vector<TH1D*> v_Gen_Top_histname1D;
                for(size_t iHist = 0; iHist < v_Gen_Top_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Gen_Top_histname1D_.at(iHist));
                    v_Gen_Top_histname1D.push_back(histname1D);
                }
                 
                //Print Histograms
                plot_2TH1(outputFolder, v_Gen_Top_histname1D.at(0), v_Gen_Top_histname1D.at(1), name, name_leg, "TopJetsPt", 4, 2, "Top Leading Jet", "Top nLeading Jet");
                plot_2TH1(outputFolder, v_Gen_Top_histname1D.at(2), v_Gen_Top_histname1D.at(3), name, name_leg, "TopJetsEta", 4, 2, "Top Leading Jet", "Top nLeading Jet");
                plot_TH1 (outputFolder, v_Gen_Top_histname1D.at(4), name, name_leg, "DeltaEta_lead_nlead_TopJets", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Gen_Top_histname1D.at(5), name, name_leg, "DeltaPhi_lead_nlead_TopJets", 4, "Top Jets", false);
                plot_TH1 (outputFolder, v_Gen_Top_histname1D.at(6), name, name_leg, "DeltaR_lead_nlead_TopJets", 4, "Top Jets", false);
                
                //---------------------                Plots for the Higgs Jets         --------------------------------------//
                //set the names for the histograms
                std::vector<TString> v_Gen_Higgs_histname1D_;{
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"leading_Higgsjet_Pt"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"subleading_Higgsjet_Pt"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"leading_Higgsjet_eta"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"subleading_Higgsjet_eta"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"DeltaEta_Higgsjet"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"DeltaPhi_Higgsjet"+ending_histname);
                    v_Gen_Higgs_histname1D_.push_back(prefix_histname+"DeltaR_Higgsjet"+ending_histname);
                }
                
                // set up histograms
                std::vector<TH1D*> v_Gen_Higgs_histname1D;
                for(size_t iHist = 0; iHist < v_Gen_Higgs_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Gen_Higgs_histname1D_.at(iHist));
                    v_Gen_Higgs_histname1D.push_back(histname1D);
                }
                   
                //Print Histograms
                plot_2TH1(outputFolder, v_Gen_Higgs_histname1D.at(0), v_Gen_Higgs_histname1D.at(1), name, name_leg, "HiggsJetsPt", 7, 6, "Higgs Leading Jet", "Higgs nLeading Jet");
                plot_2TH1(outputFolder, v_Gen_Higgs_histname1D.at(2), v_Gen_Higgs_histname1D.at(3), name, name_leg, "HiggsJetsEta", 7, 6, "Higgs Leading Jet", "Higgs nLeading Jet");
                plot_TH1 (outputFolder, v_Gen_Higgs_histname1D.at(4), name, name_leg, "DeltaEta_lead_nlead_HiggsJets", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Gen_Higgs_histname1D.at(5), name, name_leg, "DeltaPhi_lead_nlead_HiggsJets", 7, "Higgs Jets", false);
                plot_TH1 (outputFolder, v_Gen_Higgs_histname1D.at(6), name, name_leg, "DeltaR_lead_nlead_HiggsJets", 7, "Higgs Jets", false);
                   
                //--------------------       Plots for Top + Higgs Jets  ----------------------------------------//
                //set the names for the histograms
                std::vector<TString> v_Gen_TopHiggs_histname1D_;{
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Topleading_TopHiggsJets_Pt"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Topsubleading_TopHiggsJets_Pt"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Higgsleading_TopHiggsJets_Pt"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Higgssubleading_TopHiggsJets_Pt"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Topleading_TopHiggsJets_eta"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Topsubleading_TopHiggsJets_eta"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Higgsleading_TopHiggsJets_eta"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"Higgssubleading_TopHiggsJets_eta"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"MinDeltaEta"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"MinDeltaPhi"+ending_histname);
                    v_Gen_TopHiggs_histname1D_.push_back(prefix_histname+"MinDeltaR"+ending_histname);
                }
                
                // set up histograms
                std::vector<TH1D*> v_Gen_TopHiggs_histname1D;
                for(size_t iHist = 0; iHist < v_Gen_TopHiggs_histname1D_.size(); ++iHist){
                    TH1D *histname1D  = (TH1D*) file->Get(v_Gen_TopHiggs_histname1D_.at(iHist));
                    v_Gen_TopHiggs_histname1D.push_back(histname1D);
                }   
                
                //Print Histograms
                plot_4TH1(outputFolder, v_Gen_TopHiggs_histname1D.at(0), v_Gen_TopHiggs_histname1D.at(1), v_Gen_TopHiggs_histname1D.at(2), v_Gen_TopHiggs_histname1D.at(3), name, name_leg, "TopHiggsJetsPt", 4, 2, 7, 6,"Top Leading Jet", "Top nLeading Jet","Higgs Leading Jet", "Higgs nLeading Jet");
                plot_4TH1(outputFolder, v_Gen_TopHiggs_histname1D.at(4), v_Gen_TopHiggs_histname1D.at(5), v_Gen_TopHiggs_histname1D.at(6), v_Gen_TopHiggs_histname1D.at(7), name, name_leg, "TopHiggsJetsEta", 4, 2, 7, 6,"Top Leading Jet", "Top nLeading Jet","Higgs Leading Jet", "Higgs nLeading Jet");
                plot_TH1 (outputFolder, v_Gen_TopHiggs_histname1D.at(8), name, name_leg, "DeltaEta_TopHiggsJets", 8, "Top+Higgs Jets", false);
                plot_TH1 (outputFolder, v_Gen_TopHiggs_histname1D.at(9), name, name_leg, "DeltaPhi_TopHiggsJets", 8, "Top+Higgs Jets", false);
                plot_TH1 (outputFolder, v_Gen_TopHiggs_histname1D.at(10), name, name_leg, "DeltaR_TopHiggsJets", 8, "Top+Higgs Jets", false);
            }


        
            //===============       Gen-Reco Matching Plots        =================//
        
            //----------------   Plots for all Jets (No matter if they are Top or Higgs)   ---------------------------//
            //set the names for the histograms
            std::vector<TString> v_GenReco_All_histname1D_;{
                v_GenReco_All_histname1D_.push_back(prefix_histname+"RecoPtOverGenPt"+ending_histname);
                v_GenReco_All_histname1D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPt"+ending_histname);
                v_GenReco_All_histname1D_.push_back(prefix_histname+"MinDeltaR_RecoGen"+ending_histname);
            }
  
            std::vector<TString> v_GenReco_All_histname2D_;
            {
                v_GenReco_All_histname2D_.push_back(prefix_histname+"MinDeltaR_RecoGen_VS_PtGen_2D"+ending_histname);
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenP_VS_RecoP_2D"+ending_histname); 
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenPtVSGenPtMinusRecoPt_2D"+ending_histname); 
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenPtVSRecoPt_2D"+ending_histname); 
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenEtaVSRecoEta_2D"+ending_histname);  
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenPhiVSRecoPhi_2D"+ending_histname); 
                v_GenReco_All_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsMinDeltaR_2D"+ending_histname);   
                v_GenReco_All_histname2D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D"+ending_histname); 
                v_GenReco_All_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPt_2D"+ending_histname);
                v_GenReco_All_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenEta_2D"+ending_histname);
                v_GenReco_All_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPhi_2D"+ending_histname);
            }
        
            // set up histograms 
            std::vector<TH1D*> v_GenReco_All_histname1D;
            std::vector<TH2D*> v_GenReco_All_histname2D;
            for(size_t iHist = 0; iHist < v_GenReco_All_histname1D_.size(); ++iHist){
                TH1D *histname1D  = (TH1D*) file->Get(v_GenReco_All_histname1D_.at(iHist));
                v_GenReco_All_histname1D.push_back(histname1D);
            }
            for(size_t iHist = 0; iHist < v_GenReco_All_histname2D_.size(); ++iHist){
                TH2D *histname2D  = (TH2D*) file->Get(v_GenReco_All_histname2D_.at(iHist));
                v_GenReco_All_histname2D.push_back(histname2D);
            }
        
            //Print Histograms
            plot_TH1 (outputFolder, v_GenReco_All_histname1D.at(0), name, name_leg, "RecoPtOverGenPt_All", 8, "All Jets", true);
            plot_TH1 (outputFolder, v_GenReco_All_histname1D.at(1), name, name_leg, "GenPtMinusRecoPtOverGenPt_All", 8, "All Jets", true);
            plot_TH1 (outputFolder, v_GenReco_All_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_All", 8, "All Jets", true);
            plot_TH1 (outputFolder, v_GenReco_All_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_All_zoom", 8, "All Jets", true, true, 0.0, 1.0);
        
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_All_zoom","All Jets", false, 0.0, 0.0, true, 0.0, 1.0);
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(1), name, name_leg, "GenP_VS_RecoP_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(2), name, name_leg, "GenPtVSGenPtMinusRecoPt_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(3), name, name_leg, "GenPtVSRecoPt_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(4), name, name_leg, "GenEtaVSRecoEta_2D_ALl","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(5), name, name_leg, "GenPhiVSRecoPhi_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_All_zoom","All Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_All_zoom","All Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(8), name, name_leg, "RecoPtOverGenPtVsGenPt_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(9), name, name_leg, "RecoPtOverGenPtVsGenEta_2D_All","All Jets");
            plot_TH2 (outputFolder, v_GenReco_All_histname2D.at(10), name, name_leg, "RecoPtOverGenPtVsGenPhi_2D_All","All Jets");
        
            //-----------------------------  Plots for the Top Jets -----------------------------------//
            //set the names for the histograms
            std::vector<TString> v_GenReco_Top_histname1D_;{
                v_GenReco_Top_histname1D_.push_back(prefix_histname+"RecoPtOverGenPt_Topjet"+ending_histname);
                v_GenReco_Top_histname1D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPt_Topjet"+ending_histname);
                v_GenReco_Top_histname1D_.push_back(prefix_histname+"MinDeltaR_RecoGen_Topjet"+ending_histname);
            }
        
            std::vector<TString> v_GenReco_Top_histname2D_;{
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"MinDeltaR_RecoGen_VS_PtGen_Topjet_2D"+ending_histname);
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenP_VS_RecoP_Topjet_2D"+ending_histname); 
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenPtVSGenPtMinusRecoPt_Topjet_2D"+ending_histname); 
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenPtVSRecoPt_Topjet_2D"+ending_histname); 
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenEtaVSRecoEta_Topjet_2D"+ending_histname);  
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenPhiVSRecoPhi_Topjet_2D"+ending_histname); 
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsMinDeltaR_Topjet_2D"+ending_histname);   
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D"+ending_histname); 
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPt_Topjet_2D"+ending_histname);
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenEta_Topjet_2D"+ending_histname);
                v_GenReco_Top_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPhi_Topjet_2D"+ending_histname);
            }
        
            // set up histograms
            std::vector<TH1D*> v_GenReco_Top_histname1D;
            std::vector<TH2D*> v_GenReco_Top_histname2D;
            for(size_t iHist = 0; iHist < v_GenReco_Top_histname1D_.size(); ++iHist){
                TH1D *histname1D  = (TH1D*) file->Get(v_GenReco_Top_histname1D_.at(iHist));
                v_GenReco_Top_histname1D.push_back(histname1D);
            }
            for(size_t iHist = 0; iHist < v_GenReco_Top_histname2D_.size(); ++iHist){
                TH2D *histname2D  = (TH2D*) file->Get(v_GenReco_Top_histname2D_.at(iHist));
                v_GenReco_Top_histname2D.push_back(histname2D);
            }
        
            //Print Histograms 
            plot_TH1 (outputFolder, v_GenReco_Top_histname1D.at(0), name, name_leg, "RecoPtOverGenPt_Top", 4, "Top Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Top_histname1D.at(1), name, name_leg, "GenPtMinusRecoPtOverGenPt_Top", 4, "Top Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Top_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_Top", 4, "Top Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Top_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_Top_zoom", 4, "Top Jets", true, true, 0.0, 1.0);
        
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_Top","Top Jets", false, 0.0, 0.0, true, 0.0, 1.0);
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(1), name, name_leg, "GenP_VS_RecoP_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(2), name, name_leg, "GenPtVSGenPtMinusRecoPt_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(3), name, name_leg, "GenPtVSRecoPt_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(4), name, name_leg, "GenEtaVSRecoEta_2D_ALl","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(5), name, name_leg, "GenPhiVSRecoPhi_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_Top_zoom","Top Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_Top_zoom","Top Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(8), name, name_leg, "RecoPtOverGenPtVsGenPt_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(9), name, name_leg, "RecoPtOverGenPtVsGenEta_2D_Top","Top Jets");
            plot_TH2 (outputFolder, v_GenReco_Top_histname2D.at(10), name, name_leg, "RecoPtOverGenPtVsGenPhi_2D_Top","Top Jets");
        
        
        
            //-------------------------   Plots for the Higgs Jets  --------------------------------------//
            //set the names for the histograms
            std::vector<TString> v_GenReco_Higgs_histname1D_;{
                v_GenReco_Higgs_histname1D_.push_back(prefix_histname+"RecoPtOverGenPt_Higgsjet"+ending_histname);
                v_GenReco_Higgs_histname1D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPt_Higgsjet"+ending_histname);
                v_GenReco_Higgs_histname1D_.push_back(prefix_histname+"MinDeltaR_RecoGen_Higgsjet"+ending_histname);
            }
          
            std::vector<TString> v_GenReco_Higgs_histname2D_;{
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D"+ending_histname);
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenP_VS_RecoP_Higgsjet_2D"+ending_histname); 
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenPtVSGenPtMinusRecoPt_Higgsjet_2D"+ending_histname); 
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenPtVSRecoPt_Higgsjet_2D"+ending_histname); 
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenEtaVSRecoEta_Higgsjet_2D"+ending_histname);  
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenPhiVSRecoPhi_Higgsjet_2D"+ending_histname); 
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D"+ending_histname);   
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D"+ending_histname); 
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPt_Higgsjet_2D"+ending_histname);
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenEta_Higgsjet_2D"+ending_histname);
                v_GenReco_Higgs_histname2D_.push_back(prefix_histname+"RecoPtOverGenPtVsGenPhi_Higgsjet_2D"+ending_histname);
            }
        
            // set up histograms
            std::vector<TH1D*> v_GenReco_Higgs_histname1D;
            std::vector<TH2D*> v_GenReco_Higgs_histname2D;
            for(size_t iHist = 0; iHist < v_GenReco_Higgs_histname1D_.size(); ++iHist){
                TH1D *histname1D  = (TH1D*) file->Get(v_GenReco_Higgs_histname1D_.at(iHist));
                v_GenReco_Higgs_histname1D.push_back(histname1D);
            }
            for(size_t iHist = 0; iHist < v_GenReco_Higgs_histname2D_.size(); ++iHist){
                TH2D *histname2D  = (TH2D*) file->Get(v_GenReco_Higgs_histname2D_.at(iHist));
                v_GenReco_Higgs_histname2D.push_back(histname2D);
            }
        
            //Print Histograms 
            plot_TH1 (outputFolder, v_GenReco_Higgs_histname1D.at(0), name, name_leg, "RecoPtOverGenPt_Higgs", 7, "Higgs Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Higgs_histname1D.at(1), name, name_leg, "GenPtMinusRecoPtOverGenPt_Higgs", 7, "Higgs Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Higgs_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_Higgs", 7, "Higgs Jets", true);
            plot_TH1 (outputFolder, v_GenReco_Higgs_histname1D.at(2), name, name_leg, "MinDeltaR_RecoGen_Higgs", 7, "Higgs Jets", true, true, 0.0, 1.0);
        
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(0), name, name_leg, "MinDeltaR_RecoGen_VS_PtGen_2D_Higgs","Higgs Jets", false, 0.0, 0.0, true, 0.0, 1.0);
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(1), name, name_leg, "GenP_VS_RecoP_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(2), name, name_leg, "GenPtVSGenPtMinusRecoPt_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(3), name, name_leg, "GenPtVSRecoPt_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(4), name, name_leg, "GenEtaVSRecoEta_2D_ALl","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(5), name, name_leg, "GenPhiVSRecoPhi_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(6), name, name_leg, "RecoPtOverGenPtVsMinDeltaR_2D_Higgs_zoom","Higgs Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(7), name, name_leg, "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D_Higgs_zoom","Higgs Jets", true, 0.0, 1.0, false, 0.0, 0.0);
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(8), name, name_leg, "RecoPtOverGenPtVsGenPt_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(9), name, name_leg, "RecoPtOverGenPtVsGenEta_2D_Higgs","Higgs Jets");
            plot_TH2 (outputFolder, v_GenReco_Higgs_histname2D.at(10), name, name_leg, "RecoPtOverGenPtVsGenPhi_2D_Higgs","Higgs Jets");
        
            file->Close();
        }
    }
}



int main(int argc, char** argv){
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convertChannels(Channel::allowedChannelsPlotting)));    
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
    std::vector<Channel::Channel> v_channel= {
        Channel::mumu,
        Channel::emu,
        Channel::ee
    };
    
    if(opt_channel.isSet()) v_channel = Channel::convertChannels(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for (auto i_channel : v_channel)std::cout << Channel::convertChannel(i_channel) << " ";
    std::cout << "\n\n";
    
    
    
    TString prefix_histname  = "jetMatch_";
    TString prefix_histname_Mismatched  = prefix_histname+"Mismatched_";
         
    for(size_t iFileName = 0; iFileName < v_inputFileTtbar.size(); ++iFileName){
        for(const auto& iStep : v_step){
            TString fileName         = v_inputFileTtbar.at(iFileName);
            TString legFileName      = v_inputLegTtbar.at(iFileName);
            TString ending_histname  = iStep;

            printHistogram(fileName, legFileName , prefix_histname, ending_histname, false, v_channel, v_systematic);
            
            //If you want to print the histograms for the Gen-Reco Matching where at least 2 Gen Jets are matched to the same Reco Jets
            //uncomment the following:
            printHistogram(fileName, legFileName, prefix_histname_Mismatched, ending_histname, true, v_channel, v_systematic);     
        }
    } 

}
