

 #include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TF1.h>

#include "EffHist.h"

EffHist::EffHist(TString name):
h_DataNum_(0),
h_DataDeNum_(0),
h_McNum_(0),
h_McDeNum_(0),
h_effData_(0),
h_effMc_(0),
h_SF_(0),
title_(0)
{
    name_=name;
}

EffHist::EffHist(TString name, TString title,int nb,float R1,float R2)
{
    nbins=nb;
    r1=R1;
    r2=R2;
    h_n1 = new TH1F();
    h_n2 = new TH1F();
    h_n1->SetBins(nb,R1,R2);
    h_n2->SetBins(nb,R1,R2);
    h_eff = new TH1F();
    h_eff->SetNameTitle(name,title);
    h_eff->SetBins(nb,R1,R2);
    
}

EffHist::~EffHist()
{
//     delete h_n1;
//     delete h_n2;
//     delete h_eff;
}

void EffHist::setTitle(TString title)
{
    title_=title;
}

TString EffHist::getTitle()
{
    return title_;
}

void EffHist::addData(TH1D* h_DataNum, TH1D* h_DataDeNum)
{
    h_DataNum_=h_DataNum;
    h_DataDeNum_=h_DataDeNum;
}

void EffHist::addMc(TH1D* h_McNum, TH1D* h_McDeNum)
{
    h_McNum_=h_McNum;
    h_McDeNum_=h_McDeNum;
}

TH1D* EffHist::getRatio(TH1D* num, TH1D* denum,int errorType)
{
        if(errorType)
        {
            num->Sumw2();
            denum->Sumw2();
            num->Divide(num,denum,1,1,"B");
        }
        else
        {
            num->Divide(num,denum,1,1);
        }
        //cout << "Error 1 bin: " <<num->GetBinError(5) << " " << sqrt(num->GetBinContent(5)*(1-num->GetBinContent(5))/denum->GetBinContent(5)) <<endl;
        return num;
}

void EffHist::savePlotEffSF(TString savepath,char* dataEffString,char* allmcEffString,char* sfString)
{
        h_effData_=getRatio((TH1D*)(h_DataNum_->Clone()),(TH1D*)(h_DataDeNum_->Clone()),1);
        h_effMc_=getRatio((TH1D*)(h_McNum_->Clone()),(TH1D*)(h_McDeNum_->Clone()),1);
        h_SF_=getRatio((TH1D*)(h_effData_->Clone()),(TH1D*)(h_effMc_->Clone()),0);
    
    h_effData_->SetStats(0);
    h_effData_->SetTitle(title_);

    //h_effData_->GetYaxis()->SetRangeUser(0.5, 1.2);
    h_effData_->GetYaxis()->SetRangeUser(-0.01, 1.5);
    h_effData_->GetYaxis()->SetTitleOffset(1.1);
    h_effData_->SetMarkerStyle(20);
    h_effData_->SetLineColor(1);
    
    
    h_effMc_->SetLineColor(kRed);
    h_effMc_->SetMarkerColor(kRed);
    
    h_SF_->GetYaxis()->SetRangeUser(-0.01, 1.5);
    h_SF_->SetMarkerStyle(20);
    h_SF_->SetMarkerColor(4);
    h_SF_->SetLineColor(1);
    
    TLegend leg(.64,.74,.94,.94);
        leg.AddEntry(h_effData_, TString("eff data: ") + dataEffString);
        leg.AddEntry(h_effMc_, TString("eff MC: ") + allmcEffString);
//         TLegend leg(.24,.74,.54,.94);
        leg.AddEntry(h_SF_, TString("SF: ") + sfString);
        leg.SetFillColor(0);
        leg.SetTextSize(0.04);
    

    TCanvas canvas("canvas","canvas");
           h_effData_->Draw();
           h_effMc_->Draw("same");
           //TF1 *f1 = new TF1("f1", "[0]", h_SF_->GetXaxis()->GetXmin(),h_SF_->GetXaxis()->GetXmax());
//            TF1 *f1 = new TF1("f1", "[0]+[1]*x", 20,200);
//            gStyle->SetOptFit(1);
//            gStyle->SetOptStat("n");
//            h_SF_->Fit("f1","R");
           h_SF_->Draw("same");
           leg.Draw("same");
           if(savepath.Contains("ee"))DrawDecayChLabel("ee");
           if(savepath.Contains("emu"))DrawDecayChLabel("e#mu");
           if(savepath.Contains("mumu"))DrawDecayChLabel("#mu#mu");
    canvas.Print(savepath);
    
//     savepath.ReplaceAll(".pdf",4,".root",5);
//     TFile*f= new TFile(savepath,"RECREATE");
//         h_SF_->Write();
//     f->Close();
    
}


void EffHist::Filln1(double x)
{
    h_n1->Fill(x);
}
void EffHist::Filln2(double x)
{
    h_n2->Fill(x);
}

void EffHist::DrawDecayChLabel(TString decaychannel, double textSize)
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


// void EffHist::DrawEff()
// {
//      h_eff->Divide(h_n1,h_n2,1,1);
//  
//      for(int i=1;i<=nbins;i++)
//        {
//             if(h_n1->GetBinContent(i)>0) h_eff->SetBinError(i,sqrt(h_eff->GetBinContent(i)*(1-h_eff->GetBinContent(i))/h_n2->GetBinContent(i)));
//        }
//         
//        h_eff->SetMarkerSize(1.5);
//        h_eff->SetMarkerStyle(20);
//        //h_eff->Draw("pe");
//         h_eff->SetStats(0);
//         h_eff->Write();
//      
// //TFile resultfile("EffHist.root", "RECREATE");
// //     h_n1->Write();
// //     h_n2->Write();
// 
// //resultfile.Close();
//      
//      //cout << "#######:  " << h_n2->GetEntries() << endl;
// }
