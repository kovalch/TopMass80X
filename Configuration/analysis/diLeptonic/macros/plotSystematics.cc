/**
 * This macro produces a plot comparing several systematic uncertainties 
 * for a given variable in a given channel.
 * All systematic variations are normazlied to unit area.
 * 
 * Results will be store in:  systematicComparison/_channel_/_variable_.eps
 * 
 * To run this macro:
 *      root -b -l -q plotSystematics.cc++
 * 
*/



#include <string>
#include "iostream"
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPaveText.h>


const TString outdir = "systematicComparison/";


void DrawDecayChLabel(TString decaychannel, double textSize)
{
    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}


const TString xAxisTitle(TString variable)
{
    if(variable == "basic_lepton_pt_step7") return TString("p_{T}^{l} [GeV]");
    if(variable == "basic_jet_pt_step7")    return TString("p_{T}^{jet} [GeV]");
    if(variable == "HypjetMulti")           return TString("N_{jets}");
    if(variable == "HypBjetMulti")          return TString("N_{b jets}");
    if(variable == "HypToppT")              return TString("p_{T}^{t} [GeV]");
    if(variable == "HypTopRapidity")        return TString("y^{y}");
    if(variable == "HypTTBarpT")            return TString("p_{T}^{t#bar{t}} [GeV]");
    if(variable == "HypTTBarRapidity")      return TString("y^{t#bar{t}}");
    if(variable == "HypBBBarMass")          return TString("m^{b#bar{b}} [GeV]");
    if(variable == "HypLLBarMass")          return TString("m^{l#bar{l}} [GeV]");
    if(variable == "HypLLBarpT")            return TString("p_{T}^{l#bar{l}} [GeV]");
    if(variable == "HypToppTTTRestFrame")   return TString("p_{T}^{t} (t#bar{t})* [GeV]");
    if(variable == "HypTTBarDeltaPhi")      return TString("#Delta#phi(t,#bar{t})");
    return TString("");
}


void plotSystematics(TString channel, TString variable, const int rebin=1, const double xmin=0, const double xmax=0, const bool logY=0)
{
    TFile *f_nominal = TFile::Open("fullAnalysis_N007/selectionRoot/Nominal/"+channel+"/"+channel+"_ttbarsignalplustau.root");
    TFile *f_scaleup = TFile::Open("fullAnalysis_N007/selectionRoot/SCALE_UP/"+channel+"/"+channel+"_ttbarsignalplustau_scaleup.root");
    TFile *f_scaledo = TFile::Open("fullAnalysis_N007/selectionRoot/SCALE_DOWN/"+channel+"/"+channel+"_ttbarsignalplustau_scaledown.root");
    TFile *f_matchup = TFile::Open("fullAnalysis_N007/selectionRoot/MATCH_UP/"+channel+"/"+channel+"_ttbarsignalplustau_matchingup.root");
    TFile *f_matchdo = TFile::Open("fullAnalysis_N007/selectionRoot/MATCH_DOWN/"+channel+"/"+channel+"_ttbarsignalplustau_matchingdown.root");
    TFile *f_massup  = TFile::Open("fullAnalysis_N007/selectionRoot/MASS_UP/"+channel+"/"+channel+"_ttbarsignalplustau_massup.root");
    TFile *f_massdo  = TFile::Open("fullAnalysis_N007/selectionRoot/MASS_DOWN/"+channel+"/"+channel+"_ttbarsignalplustau_massdown.root");
    TFile *f_powheg  = TFile::Open("fullAnalysis_N007/selectionRoot/POWHEG/"+channel+"/"+channel+"_ttbarsignalplustau_powheg.root");
    TFile *f_mcatnlo = TFile::Open("fullAnalysis_N007/selectionRoot/MCATNLO/"+channel+"/"+channel+"_ttbarsignalplustau_mcatnlo.root");
    TFile *f_jesup   = TFile::Open("fullAnalysis_N007/selectionRoot/JES_UP/"+channel+"/"+channel+"_ttbarsignalplustau.root");
    TFile *f_jesdo   = TFile::Open("fullAnalysis_N007/selectionRoot/JES_DOWN/"+channel+"/"+channel+"_ttbarsignalplustau.root");
    
    
    TH1D *h_nominal = (TH1D*)f_nominal->Get(variable);   if(!h_nominal) {std::cout<<"Nominal: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_clone   = (TH1D*)h_nominal->Clone("h_clone");
    TH1D *h_scaleup = (TH1D*)f_scaleup->Get(variable);   if(!h_scaleup) {std::cout<<"Scale Up: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_scaledo = (TH1D*)f_scaledo->Get(variable);   if(!h_scaledo) {std::cout<<"Scale Down: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_matchup = (TH1D*)f_matchup->Get(variable);   if(!h_matchup) {std::cout<<"Matching Up: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_matchdo = (TH1D*)f_matchdo->Get(variable);   if(!h_matchdo) {std::cout<<"Matching Down: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_massup = (TH1D*)f_massup->Get(variable);   if(!h_massup) {std::cout<<"Mass Up: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_massdo = (TH1D*)f_massdo->Get(variable);   if(!h_massdo) {std::cout<<"Mass Down: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_powheg = (TH1D*)f_powheg->Get(variable);   if(!h_powheg) {std::cout<<"Powheg: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_mcatnlo = (TH1D*)f_mcatnlo->Get(variable);   if(!h_mcatnlo) {std::cout<<"MCANTLO: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_jesup = (TH1D*)f_jesup->Get(variable);   if(!h_jesup) {std::cout<<"JES Up: No histogram   "<<variable<<std::endl; exit(11);};
    TH1D *h_jesdo = (TH1D*)f_jesdo->Get(variable);   if(!h_jesdo) {std::cout<<"JES Down: No histogram   "<<variable<<std::endl; exit(11);};
    
    h_nominal->SetLineColor(kBlack); h_nominal->SetLineStyle(1);  h_nominal->Rebin(rebin);
    h_clone->SetLineColor(0);        h_clone->SetLineStyle(0);    h_clone->Rebin(rebin);
    h_scaleup->SetLineColor(kRed);   h_scaleup->SetLineStyle(1);  h_scaleup->Rebin(rebin);
    h_scaledo->SetLineColor(kRed);   h_scaledo->SetLineStyle(3);  h_scaledo->Rebin(rebin);
    h_matchup->SetLineColor(kBlue);  h_matchup->SetLineStyle(1);  h_matchup->Rebin(rebin);
    h_matchdo->SetLineColor(kBlue);  h_matchdo->SetLineStyle(3);  h_matchdo->Rebin(rebin);
    h_massup->SetLineColor(kGreen);  h_massup->SetLineStyle(1);   h_massup->Rebin(rebin);
    h_massdo->SetLineColor(kGreen);  h_massdo->SetLineStyle(3);   h_massdo->Rebin(rebin);
    h_powheg->SetLineColor(40);      h_powheg->SetLineStyle(1);   h_powheg->Rebin(rebin);
    h_mcatnlo->SetLineColor(40);     h_mcatnlo->SetLineStyle(3);  h_mcatnlo->Rebin(rebin);
    h_jesup->SetLineColor(12);     h_jesup->SetLineStyle(1);   h_jesup->Rebin(rebin);
    h_jesdo->SetLineColor(12);     h_jesdo->SetLineStyle(3);   h_jesdo->Rebin(rebin);
    
    h_nominal->SetTitle("");
    h_nominal->GetYaxis()->SetTitleOffset(1.4);
    h_nominal->SetMaximum(1.25*h_nominal->GetMaximum());
    h_nominal->GetXaxis()->SetTitle(xAxisTitle(variable));
    if(xmin!=0 && xmax !=0) h_nominal->GetXaxis()->SetRangeUser(xmin, xmax);
    h_nominal->GetYaxis()->SetTitle("a.u.");
    
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetFillColor(0); leg->SetFillStyle(0); 
    leg->SetNColumns(2);
    leg->AddEntry(h_nominal, "Nominal", "l");     leg->AddEntry(h_clone, " ", "l");
    leg->AddEntry(h_scaleup, "Scale Up", "l");    leg->AddEntry(h_scaledo, "Scale Down", "l");
    leg->AddEntry(h_matchup, "Matching Up", "l"); leg->AddEntry(h_matchdo, "Matching Down", "l");
    leg->AddEntry(h_massup, "Mass Up", "l");      leg->AddEntry(h_massdo, "Mass Down", "l");
    leg->AddEntry(h_jesup, "JES Up", "l");      leg->AddEntry(h_jesdo, "JES Down", "l");
    leg->AddEntry(h_powheg, "Had. Up", "l");      leg->AddEntry(h_mcatnlo, "Had. Down", "l");

    TCanvas *c = new TCanvas("", "", 0, 0, 1000, 1000);
    c->SetLogy(logY);
    h_nominal->DrawNormalized();
    h_scaleup->DrawNormalized("same");  h_scaledo->DrawNormalized("same");
    h_matchup->DrawNormalized("same");  h_matchdo->DrawNormalized("same");
    h_massup->DrawNormalized("same");   h_massdo->DrawNormalized("same");
    h_jesup->DrawNormalized("same");    h_jesdo->DrawNormalized("same");
    h_powheg->DrawNormalized("same");   h_mcatnlo->DrawNormalized("same");
    DrawDecayChLabel(channel, 0.04);
    leg->Draw("same");
    c->Print(outdir+channel+"/"+variable+".eps");
    delete c;

    
    delete h_nominal;
    delete h_scaleup; delete h_scaledo;
    delete h_massup;  delete h_massdo;
    delete h_matchup; delete h_matchdo;
    delete h_jesup;   delete h_jesdo;
    delete h_powheg;  delete h_mcatnlo;
    
    f_nominal->Close();
    f_scaleup->Close(); f_scaledo->Close();
    f_matchup->Close(); f_matchdo->Close();
    f_jesup->Close();   f_jesdo->Close();
    f_massup->Close();  f_massdo->Close();
    f_powheg->Close();  f_mcatnlo->Close();
}







void plotSystematics()
{
    const std::string channels[3] = {"ee", "emu", "mumu"};
    
    gStyle->SetOptStat(0);

    for(size_t iter = 0 ; iter<3; iter++){
        gSystem->mkdir(outdir+channels[iter], kTRUE);
        plotSystematics(channels[iter], "basic_lepton_pt_step7", 2);
        plotSystematics(channels[iter], "basic_jet_pt_step7", 2);
        plotSystematics(channels[iter], "HypjetMulti", 1, 1.5, 8.4);
        plotSystematics(channels[iter], "HypBjetMulti", 1, 0.5, 6.4);
        plotSystematics(channels[iter], "HypToppT", 20);
        plotSystematics(channels[iter], "HypTopRapidity", 5, -2.4, 2.4);
        plotSystematics(channels[iter], "HypTTBarpT", 20);
        plotSystematics(channels[iter], "HypTTBarRapidity", 5, -2.4, 2.4);
//         plotSystematics(channels[iter], "HypBBBarMass");
//         plotSystematics(channels[iter], "HypLLBarMass");
//         plotSystematics(channels[iter], "HypLLBarpT");
//         plotSystematics(channels[iter], "HypToppTTTRestFrame");
//         plotSystematics(channels[iter], "HypTTBarDeltaPhi");
    }
}


