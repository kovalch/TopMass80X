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
#include "THStack.h"
#include "TPad.h"
#include "TLegend.h"

#include "ProgramOptionsReader.h"
#include "Helper.h"

#include <iostream>
#include <boost/range/adaptor/reversed.hpp>

typedef ProgramOptionsReader po;

TopMassControlPlots::TopMassControlPlots()
{
  doPlots();
}

void TopMassControlPlots::doPlots()
{
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  Helper* helper = new Helper();
  helper->SetTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
	gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(1.5);

  std::string path(po::GetOption<std::string>("analysisConfig.samplePath"));

  // DUMMY to get TTree structure
  TChain* myChain = new TChain("analyzeKinFit/eventTree");
  myChain->Add((path+std::string("Z2_S12_P11TeV_sig.root")).c_str());

  hists.push_back(MyHistogram("1", "top.fitTop1[0].M()", myChain, ";m_{t}^{fit} [GeV]; Events", 75, 50, 400));
  /*
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
  */
  
  samples.push_back(MySample("Data", "MJP12*v1_data.root", kData, kBlack));
  samples.push_back(MySample("t#bar{t}", "Z2_S12_ABS_JES_100_169_5_sig.root", kSig, kRed+1));
  samples.push_back(MySample("QCD", "QCDMixing_MJPS12*v1_data.root", kBkg, kYellow));

  TTreeFormula weight("weight", po::GetOption<std::string>("weight").c_str(), myChain);
  TTreeFormula sel   ("sel"   , po::GetOption<std::string>("analysisConfig.selection").c_str(), myChain);
  TTreeFormula selCP ("selCP" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" && ")+po::GetOption<std::string>("analysisConfig.selectionCP")).c_str(), myChain);
  TTreeFormula selWP ("selWP" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" && ")+po::GetOption<std::string>("analysisConfig.selectionWP")).c_str(), myChain);
  TTreeFormula selUN ("selUN" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" && ")+po::GetOption<std::string>("analysisConfig.selectionUN")).c_str(), myChain);

  delete myChain;

  {
    int bkgCounter = -1;
    for(MySample& sample : samples){
      if(sample.type == kBkg ) ++bkgCounter;
      TChain* chain = new TChain("analyzeKinFit/eventTree");
      int nFiles = chain->Add((path+sample.file).c_str());
      std::cout << "Adding " << nFiles << " files for " << sample.name << " (" << sample.file << ")" << std::endl;

      for(MyHistogram& hist : hists) {
        hist.var->SetTree(chain);
        if     (sample.type == kSig ) hist.AddSignal(&sample);
        else if(sample.type == kBkg ) hist.AddBackground(&sample);
      }
      weight.SetTree(chain);
      sel   .SetTree(chain);
      selCP .SetTree(chain);
      selWP .SetTree(chain);
      selUN .SetTree(chain);

      for(int i = 0; ; ++i){
        long entry = chain->LoadTree(i);
        if(entry  < 0) break;
        if(entry == 0){
          for(auto& hist : hists)
            hist.var->UpdateFormulaLeaves();
          weight.UpdateFormulaLeaves();
          sel   .UpdateFormulaLeaves();
          selCP .UpdateFormulaLeaves();
          selWP .UpdateFormulaLeaves();
          selUN .UpdateFormulaLeaves();
        }
        if(!weight.GetNdata()) continue;
        if(!sel   .GetNdata()) continue;

        for(int j = 0, l = sel.GetNdata(); j < l; ++j){
          if(!sel.EvalInstance(j)) continue;
          for(MyHistogram& hist : hists){
            if(!hist.var->GetNdata()) continue;
            if     (sample.type == kData) hist.data           ->Fill(hist.var->EvalInstance(j));
            else if(sample.type == kSig ) {
              if(selCP.EvalInstance(j)) hist.sig[2]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j));
              if(selWP.EvalInstance(j)) hist.sig[1]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j));
              if(selUN.EvalInstance(j)) hist.sig[0]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j));
            }
            else if(sample.type == kBkg ) hist.bkg[bkgCounter]->Fill(hist.var->EvalInstance(j));
          }
        }
      }
    }
  }
  for(MyHistogram& hist : hists){
    int bins = hist.data->GetNbinsX()+1;
    double integralD = hist.data->Integral(0,bins);
    double integralS = 0; for(TH1F* sig : hist.sig) integralS += sig->Integral(0,bins);
    double integralB = 0; for(TH1F* bkg : hist.bkg) integralB += bkg->Integral(0,bins);
    double fSig = po::GetOption<double>("templates.fSig");
    for(TH1F* sig : hist.sig) sig->Scale(    fSig *integralD/integralS);
    for(TH1F* bkg : hist.bkg) bkg->Scale((1.-fSig)*integralD/integralB);
    TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
    canv->cd();
    
    hist.data->GetYaxis()->SetRangeUser(1, hist.data->GetMaximum()*1.6);
    hist.data->Draw("p");
    THStack* stack = new THStack("stack", "");
    for(TH1F* bkg : boost::adaptors::reverse(hist.bkg)) stack->Add(bkg);
    for(TH1F* sig : boost::adaptors::reverse(hist.sig)) stack->Add(sig);
    stack    ->Draw("hist same");
    hist.data->Draw("p same");
    
    TLegend* leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
    leg0->SetTextSize(0.03);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->AddEntry( hist.data, hist.data->GetTitle(), "P" );
    for(TH1F* sig : hist.sig) leg0->AddEntry( sig, sig->GetTitle(), "F" );
    for(TH1F* bkg : hist.bkg) leg0->AddEntry( bkg, bkg->GetTitle(), "F" );
    leg0->Draw();
    
    gPad->RedrawAxis();
    
    canv->Print((std::string("plot/controlplots/")+std::string(hist.data->GetName())+std::string(".eps")).c_str(),"eps");
    canv->Print((std::string("plot/controlplots/")+std::string(hist.data->GetName())+std::string(".png")).c_str(),"png");
  }
}
