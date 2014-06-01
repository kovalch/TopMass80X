#ifndef EffHist_h
#define EffHist_h



#include <TH1.h>
#include <TH1F.h>
#include <iostream>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>
#include <TStyle.h>

using namespace std;




class EffHist {

    public:
       EffHist(TString name);
       EffHist(TString name,TString title,int nb,float R1,float R2);
       ~EffHist();
       void Filln1(double x);
       void Filln2(double x);
       void DrawEff();
       void addMc(TH1D* h_McNum,TH1D* h_McDeNum);
       void addData(TH1D* h_DataNum,TH1D* h_DataDeNum);
       TH1D *getRatio(TH1D *num, TH1D *demun, int errorType);
       void savePlotEffSF(TString savepath,char* dataEffString,char* allmcEffString,char* sfString);
       void DrawDecayChLabel(TString decaychannel, double textSize= 0.04);
       void setTitle(TString title);
       TString getTitle();
       
    private:
       int nbins;
       float r1;
       float r2;
       TH1F *h_n1;
       TH1F *h_n2;
       TH1F *h_eff;
       
       TH1D *h_DataNum_;
       TH1D *h_DataDeNum_;
       TH1D *h_McNum_;
       TH1D *h_McDeNum_;
       
       TH1D *h_effData_;
       TH1D *h_effMc_;
       TH1D *h_SF_;
       
       TString name_;
       TString title_;

        
  };


#endif
