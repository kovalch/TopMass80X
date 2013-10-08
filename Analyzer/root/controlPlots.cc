#include "controlPlots.h"

#include "TCanvas.h"
#include "TDirectory.h"
//#include "TLegend.h"
//#include "TROOT.h"
//#include "TString.h"
#include "TStyle.h"
//#include "TTree.h"
//#include "TLorentzVector.h"
//#include "TH2F.h"
//#include "TSeqCollection.h"
//#include "TGraphAsymmErrors.h"
//#include "TPaveLabel.h"

#include "ProgramOptionsReader.h"

#include <iostream>

typedef ProgramOptionsReader po;

TopMassControlPlots::TopMassControlPlots()
{
  doPlots();
}

void TopMassControlPlots::doPlots()
{
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  std::string path(po::GetOption<std::string>("analysisConfig.samplePath"));

  // DUMMY to get TTree structure
  TChain* myChain = new TChain("analyzeKinFit/eventTree");
  myChain->Add((path+std::string("Z2_S12_P11TeV_sig.root")).c_str());

  hists.push_back(MyHistogram("1", "top.fitTop1[0].M()", myChain, ";m_{t}^{fit} [GeV]; Events", 125, 100, 350));
  hists.push_back(MyHistogram("2", "(top.recoW1[0].M()+top.recoW2[0].M())/2.0", myChain, ";m_{W}^{reco} [GeV]; Events", 60, 65, 125));
  hists.push_back(MyHistogram("3", "top.fitProb[0]", myChain, ";P_{gof}; Events", 50, 0, 1.0));
  hists.push_back(MyHistogram("4", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", myChain, ";#DeltaR_{b#bar{b}}; Events", 50, 1, 6));
  hists.push_back(MyHistogram("5", "jet.jet[0].Pt()", myChain, ";p_{T}^{1} [GeV]; Events", 45, 0, 450));
  hists.push_back(MyHistogram("6", "jet.jet[1].Pt()", myChain, ";p_{T}^{2} [GeV]; Events", 40, 0, 400));
  hists.push_back(MyHistogram("7", "jet.jet[2].Pt()", myChain, ";p_{T}^{3} [GeV]; Events", 50, 0, 250));
  hists.push_back(MyHistogram("8", "jet.jet[3].Pt()", myChain, ";p_{T}^{4} [GeV]; Events", 40, 0, 200));
  hists.push_back(MyHistogram("9", "jet.jet[4].Pt()", myChain, ";p_{T}^{5} [GeV]; Events", 50, 0, 150));
  hists.push_back(MyHistogram("0", "jet.jet[5].Pt()", myChain, ";p_{T}^{6} [GeV]; Events", 50, 0, 100));
  hists.push_back(MyHistogram("A", "top.fitTop1[0].Pt()", myChain, ";p_{T,t}^{fit} [GeV]; Events", 100, 0, 1000));

  TTreeFormula weight("weight", po::GetOption<std::string>("weight").c_str(), myChain);
  TTreeFormula sel   ("sel"   , po::GetOption<std::string>("analysisConfig.selection").c_str(), myChain);

  delete myChain;

  for(MyHistogram& hist : hists){
    hist.AddSignal();
    hist.AddBackground();
  }
  {
    int counter = -1;
    for(std::string file : {"MJP12*v1_data.root", "Z2_S12_ABS_JES_100_172_5_sig.root", "QCDMixing_MJPS12*v1_data.root"}){
      ++counter;

      TChain* chain = new TChain("analyzeKinFit/eventTree");
      int nFiles = chain->Add((path+file).c_str());
      std::cout << "Adding " << nFiles << " files for " << file << std::endl;

      for(MyHistogram& hist : hists)
        hist.var->SetTree(chain);
      weight.SetTree(chain);
      sel   .SetTree(chain);

      for(int i = 0; ; ++i){
        long entry = chain->LoadTree(i);
        if(entry  < 0) break;
        if(entry == 0){
          for(auto& hist : hists)
            hist.var->UpdateFormulaLeaves();
          weight.UpdateFormulaLeaves();
          sel   .UpdateFormulaLeaves();
        }
        if(!weight.GetNdata()) continue;
        if(!sel   .GetNdata()) continue;
        for(int j = 0, l = sel.GetNdata(); j < l; ++j){
          if(!sel.EvalInstance(j)) continue;
          for(MyHistogram& hist : hists){
            if(!hist.var->GetNdata()) continue;
            if     (counter == 0) hist.data  ->Fill(hist.var->EvalInstance(j));
            else if(counter == 1) hist.sig[0]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j));
            else if(counter == 2) hist.bkg[0]->Fill(hist.var->EvalInstance(j));
          }
        }
      }
    }
  }
  for(MyHistogram& hist : hists){
    int bins = hist.data->GetNbinsX()+1;
    double integralD = hist.data->Integral(0,bins);
    double integralS = hist.sig[0]->Integral(0,bins);
    double integralB = hist.bkg[0]->Integral(0,bins);
    double fSig = po::GetOption<double>("templates.fSig");
    hist.sig[0]->Scale(    fSig *integralD/integralS);
    hist.bkg[0]->Scale((1.-fSig)*integralD/integralB);
    TCanvas* canv = new TCanvas();
    canv->cd();
    hist.data  ->Draw("p");
    hist.sig[0]->Draw("hist same");
    hist.bkg[0]->Draw("hist same");
    TH1F* comb = (TH1F*)hist.sig[0]->Clone();
    comb->Add(hist.bkg[0]);
    comb->SetLineColor(8);
    comb->SetLineStyle(0);
    comb->Draw("hist same");
    canv->Print((std::string(hist.data->GetName())+std::string(".eps")).c_str(),"eps");
    canv->Print((std::string(hist.data->GetName())+std::string(".png")).c_str(),"png");
  }
}
