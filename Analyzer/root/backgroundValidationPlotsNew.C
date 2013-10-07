#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
//#include "TLegend.h"
//#include "TROOT.h"
//#include "TString.h"
#include "TStyle.h"
//#include "TTree.h"
#include "TTreeFormula.h"
//#include "TLorentzVector.h"
//#include "TClonesArray.h"
#include "TH1F.h"
//#include "TH2F.h"
//#include "TSeqCollection.h"
//#include "TGraphAsymmErrors.h"
//#include "TPaveLabel.h"

#include <iostream>

struct MyHistogram{

  MyHistogram(std::string name, std::string formula, TChain* chain, std::string title, int nBins, double min, double max) :
    var(new TTreeFormula((std::string("f")+name).c_str(), formula.c_str(), chain)),
    data(new TH1F((std::string("hD")+name).c_str(), title.c_str(), nBins, min, max)),
    sig(new TH1F((std::string("hS")+name).c_str(), title.c_str(), nBins, min, max)),
    bkg(new TH1F((std::string("hB")+name).c_str(), title.c_str(), nBins, min, max))
  {
    data->SetLineWidth(2);
    data->SetLineColor(kBlack);
    data->SetMarkerStyle(20);
    data->SetMarkerColor(kBlack);

    sig->SetLineWidth(2);
    sig->SetLineColor(kRed);
    sig->SetLineStyle(9);

    bkg->SetLineWidth(2);
    bkg->SetLineColor(kBlue);
    bkg->SetLineStyle(2);
  }

  TTreeFormula* var;
  TH1F *data, *sig, *bkg;
};

void backgroundValidationPlotsNew()
{
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  std::string path("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/2012/02/");

  // DUMMY to get TTree structure
  TChain* myChain = new TChain("analyzeKinFit/eventTree");
  myChain->Add((path+std::string("Z2_S12_Lep_P11TeV_sig.root")).c_str());

  std::vector<MyHistogram> hists;
  hists.push_back(MyHistogram("1", "top.fitTop1[0].M()", myChain, ";m_{t}^{fit} [GeV]; Events", 125, 100, 350));
  hists.push_back(MyHistogram("2", "(top.recoW1[0].M()+top.recoW2[0].M())/2.0", myChain, ";m_{W}^{reco} [GeV]; Events", 60, 65, 125));
  hists.push_back(MyHistogram("3", "top.fitProb[0]", myChain, ";P_{gof}; Events", 50, 0, 1.0));
  hists.push_back(MyHistogram("4", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", myChain, ";#DeltaR_{b#bar{b}}; Events", 50, 1, 6));
  hists.push_back(MyHistogram("5", "jet.jet[0].Pt()", myChain, ";p_{T}^{1}; Events", 45, 0, 450));
  hists.push_back(MyHistogram("6", "jet.jet[1].Pt()", myChain, ";p_{T}^{2}; Events", 40, 0, 400));
  hists.push_back(MyHistogram("7", "jet.jet[2].Pt()", myChain, ";p_{T}^{3}; Events", 50, 0, 250));
  hists.push_back(MyHistogram("8", "jet.jet[3].Pt()", myChain, ";p_{T}^{4}; Events", 40, 0, 200));
  hists.push_back(MyHistogram("9", "jet.jet[4].Pt()", myChain, ";p_{T}^{5}; Events", 50, 0, 150));
  hists.push_back(MyHistogram("0", "jet.jet[5].Pt()", myChain, ";p_{T}^{6}; Events", 50, 0, 100));
  hists.push_back(MyHistogram("A", "top.fitTop1[0].Pt()", myChain, ";p_{T,t}^{fit} [GeV]; Events", 100, 0, 1000));

  TTreeFormula weight = TTreeFormula("weight", "weight.combinedWeight", myChain);
  TTreeFormula sel    = TTreeFormula("sel"   , "top.fitProb[0] > 0.1 && (pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2)) > 1.5*1.5", myChain);
  TTreeFormula index  = TTreeFormula("index" , "Iteration$", myChain);

  delete myChain;

  {
    int counter = -1;
    for(std::string file : {"MJP12*_v1_data.root", "Z2_S12_Had*_ABS_JES_100_172_5_sig.root", "QCDMixing_MJPS12*_v1_data.root"}){
      ++counter;

      TChain* chain = new TChain("analyzeKinFit/eventTree");
      int nFiles = chain->Add((path+file).c_str());
      std::cout << "Adding " << nFiles << " files for " << file << std::endl;

      for(auto& hist : hists)
        hist.var->SetTree(chain);
      weight.SetTree(chain);
      sel   .SetTree(chain);
      index .SetTree(chain);

      for(int i = 0; ; ++i){
        long entry = chain->LoadTree(i);
        if(entry  < 0) break;
        if(entry == 0){
          for(auto& hist : hists)
            hist.var->UpdateFormulaLeaves();
          weight.UpdateFormulaLeaves();
          sel   .UpdateFormulaLeaves();
          index .UpdateFormulaLeaves();
        }
        if(!sel.GetNdata()) continue;
        if(!sel.EvalInstance(0)) continue;
        for(int j = 0, l = index.GetNdata(); j < l; ++j){
          if(!sel.EvalInstance(j)) continue;
          for(auto& hist : hists){
            if     (counter == 0) hist.data->Fill(hist.var->EvalInstance(j));
            else if(counter == 1) hist.sig ->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j));
            else if(counter == 2) hist.bkg ->Fill(hist.var->EvalInstance(j));
          }
        }
      }
    }
  }
  for(auto& hist : hists){
    int bins = hist.data->GetNbinsX()+1;
    double integralD = hist.data->Integral(0,bins);
    double integralS = hist.sig ->Integral(0,bins);
    double integralB = hist.bkg ->Integral(0,bins);
    double fSig = 0.460272275;
    hist.sig->Scale(    fSig *integralD/integralS);
    hist.bkg->Scale((1.-fSig)*integralD/integralB);
    TCanvas* canv = new TCanvas();
    canv->cd();
    hist.data->Draw("p");
    hist.sig ->Draw("hist same");
    hist.bkg ->Draw("hist same");
    TH1F* comb = (TH1F*)hist.sig->Clone();
    comb->Add(hist.bkg);
    comb->SetLineColor(8);
    comb->SetLineStyle(0);
    comb->Draw("hist same");
  }
}
