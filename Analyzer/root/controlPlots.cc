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
  std::string channel(po::GetOption<std::string>("channel"));
  int channelID_ = Helper::channelID();
  double lumi    = po::GetOption<double>("lumi");

  // DUMMY to get TTree structure
  TChain* myChain = new TChain("analyzeKinFit/eventTree");
  myChain->Add("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/2012/Skim_02/Z2_S12_P11TeV_sig.root");

  // masses
  hists.push_back(MyHistogram("fitTop1Mass", "top.fitTop1.M()", myChain, ";m_{t}^{fit} [GeV]; Permutations", 75, 50, 400));
  hists.push_back(MyHistogram("recoWAveMass", "(top.recoW1[0].M()+top.recoW2[0].M())/2.0", myChain, ";m_{W}^{reco} [GeV]; Events", 60, 65, 125));
  hists.push_back(MyHistogram("recoW1Mass", "top.recoW1.M()", myChain, ";m_{W}^{reco} [GeV]; Permutations", 60, 0, 300));
  
  // light pulls
  hists.push_back(MyHistogram("pullW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pull[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", myChain, ";#theta_{pull}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", myChain, ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", myChain, ";#theta_{pull,ch}^{q_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod2_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),-(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta())))))", myChain, ";#theta_{pull,ch}^{q_{2} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", myChain, ";#theta_{pull,ch}^{q_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  
  // b pulls
  hists.push_back(MyHistogram("pullChargedB1_B2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))", myChain, ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", myChain, ";#theta_{pull,ch}^{b_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB2_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitB2.Phi()),-(top.fitB1.Eta()-top.fitB2.Eta())))))", myChain, ";#theta_{pull,ch}^{b_{2} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", myChain, ";#theta_{pull,ch}^{b_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  
  // mixed pulls
  hists.push_back(MyHistogram("pullChargedW1Prod1_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))", myChain, ";#theta_{pull}^{q_{1} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitB1.Phi()),-(top.fitW1Prod1.Eta()-top.fitB1.Eta())))))", myChain, ";#theta_{pull}^{b_{1} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));
  
  // others
  hists.push_back(MyHistogram("fitProb", "top.fitProb", myChain, ";P_{gof}; Permutations", 50, 0, 1.0));
  hists.push_back(MyHistogram("fitProbBest", "top.fitProb[0]", myChain, ";P_{gof}; Events", 50, 0, 1.0));
  hists.push_back(MyHistogram("deltaRbbBest", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", myChain, ";#DeltaR_{b#bar{b}}; Events", 50, 1, 6));
  
  // pts
  hists.push_back(MyHistogram("jet1Pt", "jet.jet[0].Pt()", myChain, ";p_{T}^{1} [GeV]; Events", 45, 0, 450));
  hists.push_back(MyHistogram("jet2Pt", "jet.jet[1].Pt()", myChain, ";p_{T}^{2} [GeV]; Events", 40, 0, 400));
  hists.push_back(MyHistogram("jet3Pt", "jet.jet[2].Pt()", myChain, ";p_{T}^{3} [GeV]; Events", 50, 0, 250));
  hists.push_back(MyHistogram("jet4Pt", "jet.jet[3].Pt()", myChain, ";p_{T}^{4} [GeV]; Events", 40, 0, 200));
  hists.push_back(MyHistogram("jet5Pt", "jet.jet[4].Pt()", myChain, ";p_{T}^{5} [GeV]; Events", 50, 0, 150));
  hists.push_back(MyHistogram("jet6Pt", "jet.jet[5].Pt()", myChain, ";p_{T}^{6} [GeV]; Events", 50, 0, 100));
  hists.push_back(MyHistogram("topPt", "top.fitTop1.Pt()", myChain, ";p_{T,t}^{fit} [GeV]; Permutations", 40, 0, 400));
  hists.push_back(MyHistogram("fitRlb", "(top.fitB1.Pt()+top.fitB2.Pt())/(top.fitW1Prod1.Pt()+top.fitW1Prod2.Pt())", myChain, ";R_{lb}^{fit}; Permutations", 40, 0, 4));
  hists.push_back(MyHistogram("recoRlb", "(top.recoB1.Pt()+top.recoB2.Pt())/(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())", myChain, ";R_{lb}^{reco}; Permutations", 40, 0, 4));
  
  // datasets
  if (channelID_ == Helper::kAllJets) {
    samples.push_back(MySample("Data", "MJP12*v1_data", kData, kBlack));
    samples.push_back(MySample("t#bar{t}", "Z2_S12_ABS_JES_100_172_5_sig", kSig, kRed+1));
    samples.push_back(MySample("QCD", "QCDMixing_MJPS12*v1_data", kBkg, kYellow));
  }
  else {
    samples.push_back(MySample("Data", "Run2012", kData, kBlack));
    samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kSig, kRed+1, 1, lumi/1000.));
    samples.push_back(MySample("Z+Jets", "Summer12_ZJets", kBkg, kAzure-2, 1, lumi/1000.));
    samples.push_back(MySample("W+Jets", "Summer12_WJets", kBkg, kGreen-3, 1, lumi/1000.));
    samples.push_back(MySample("single top", "Summer12_singleTop", kBkg, kMagenta, 1, lumi/1000.));
    //* Signal variations
    samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_MGDecays", kSigVar, kRed+1, 1, lumi/1000.*9./4.));
    samples.push_back(MySample("t#bar{t}, P11", "Summer12_TTJets1725_MGDecays_P11", kSigVar, kMagenta+1, 9, lumi/1000.*9./4.));
    samples.push_back(MySample("t#bar{t}, P11noCR", "Summer12_TTJets1725_MGDecays_P11noCR", kSigVar, kCyan+1, 2, lumi/1000.*9./4.));//*/
  }
  

  TTreeFormula weight("weight", po::GetOption<std::string>("weight").c_str(), myChain);
  TTreeFormula sel   ("sel"   , po::GetOption<std::string>("analysisConfig.selection").c_str(), myChain);
  TTreeFormula selCP ("selCP" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionCP")).c_str(), myChain);
  TTreeFormula selWP ("selWP" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionWP")).c_str(), myChain);
  TTreeFormula selUN ("selUN" , (po::GetOption<std::string>("analysisConfig.selection")
            +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionUN")).c_str(), myChain);

  delete myChain;

  {
    int bkgCounter = -1;
    int sigVarCounter = -1;
    for(MySample& sample : samples){
      if(sample.type == kBkg    ) ++bkgCounter;
      if(sample.type == kSigVar ) ++sigVarCounter;
      TChain* chain; int nFiles = 0;
      if (channelID_ == Helper::kAllJets) {
        chain = new TChain("analyzeKinFit/eventTree");
        nFiles = chain->Add((path+sample.file+std::string(".root")).c_str());
      }
      else {
        chain = new TChain("analyzeHitFit/eventTree");
        if (channelID_ != Helper::kElectronJets) nFiles += chain->Add((path+sample.file+std::string("_muon/job_*.root")).c_str());
        if (channelID_ != Helper::kMuonJets) nFiles += chain->Add((path+sample.file+std::string("_electron/job_*.root")).c_str());
      }
      std::cout << "Adding " << nFiles << " files for " << sample.name << " (" << sample.file << ")" << std::endl;

      for(MyHistogram& hist : hists) {
        hist.Init(chain);
        if     (sample.type == kSig    ) hist.AddSignal(&sample);
        else if(sample.type == kBkg    ) hist.AddBackground(&sample);
        else if(sample.type == kSigVar ) hist.AddSignalVariation(&sample);
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
        selCP.GetNdata(); selWP.GetNdata(); selUN.GetNdata();

        for(int j = 0, l = sel.GetNdata(); j < l; ++j){
          if(!sel.EvalInstance(j)) continue;
          for(MyHistogram& hist : hists){
            if(!hist.var->GetNdata()) continue;
            if     (sample.type == kData) {
              hist.data->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kSig ) {
              if     (selCP.EvalInstance(j)) hist.sig[2]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selWP.EvalInstance(j)) hist.sig[1]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selUN.EvalInstance(j)) hist.sig[0]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kBkg ) {
              hist.bkg[bkgCounter]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              // TODO: Background normalization in signal variations not correct for AllJets
              for(TH1F* sigvar : hist.sigvar) sigvar->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kSigVar ) hist.sigvar[sigVarCounter]->Fill(hist.var->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
          }
        }
      }
      for(MyHistogram& hist : hists) {
        hist.var->Clear();
      }
    }
  }
  bool firstHist = true;
  for(MyHistogram& hist : hists){
    int bins = hist.data->GetNbinsX()+1;
    if (firstHist) std::cout << "Yields" << std::endl;
    double integralD = hist.data->Integral(0,bins);
    if (firstHist) std::cout << "Data:       " << integralD << std::endl;
    double integralS = 0; for(TH1F* sig : hist.sig) {
      if (firstHist) std::cout << "  " << sig->GetTitle() << ": " << sig->Integral(0,bins) << std::endl;
      integralS += sig->Integral(0,bins);
    }
    if (firstHist) std::cout << "Signal:     " << integralS << std::endl;
    double integralB = 0; for(TH1F* bkg : hist.bkg) {
      if (firstHist) std::cout << "  " << bkg->GetTitle() << ": " << bkg->Integral(0,bins) << std::endl;
      integralB += bkg->Integral(0,bins);
    }
    if (firstHist) std::cout << "Background: " << integralB << std::endl;
    if (firstHist) std::cout << "MC Total:   " << integralS+integralB << std::endl;
    firstHist = false;
    
    if (channelID_ == Helper::kAllJets) {
      double fSig = po::GetOption<double>("templates.fSig");
      for(TH1F* sig    : hist.sig)    sig->Scale(    fSig *integralD/integralS);
      for(TH1F* bkg    : hist.bkg)    bkg->Scale((1.-fSig)*integralD/integralB);
      for(TH1F* sigvar : hist.sigvar) sigvar->Scale(integralD/sigvar->Integral(0,bins));
    }
    else {
      for(TH1F* sig    : hist.sig)    sig->Scale(integralD/(integralS+integralB));
      for(TH1F* bkg    : hist.bkg)    bkg->Scale(integralD/(integralS+integralB));
      for(TH1F* sigvar : hist.sigvar) sigvar->Scale(integralD/sigvar->Integral(0,bins));
      
      // Double bin error due to correlations
      for (int i = 0; i < hist.data->GetNbinsX(); i++) {
        hist.data->SetBinError(i, hist.data->GetBinError(i) * 2);
      }
    }
    
    // Draw control plot
    
    TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
    canv->cd();
    
    hist.data->GetYaxis()->SetRangeUser(1, hist.data->GetMaximum()*1.5);
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
    for(TH1F* sig : hist.sig) leg0->AddEntry( sig, sig->GetTitle(), "F" );
    leg0->Draw();
    
    TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    for(TH1F* bkg : hist.bkg) leg1->AddEntry( bkg, bkg->GetTitle(), "F" );
    leg1->AddEntry( hist.data, hist.data->GetTitle(), "P" );
    leg1->Draw();
    
    helper->DrawCMS();
    
    gPad->RedrawAxis();
    
    canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.data->GetName())+std::string(".eps")).c_str(),"eps");
    canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.data->GetName())+std::string(".png")).c_str(),"png");
    
    // Draw signal variation plot
    
    hist.data->GetYaxis()->UnZoom();
    hist.data->GetYaxis()->SetRangeUser(hist.data->GetMinimum()*0.9, hist.data->GetMaximum()*1.15);
    hist.data->Draw("p");
    for(TH1F* sigvar : hist.sigvar) sigvar->Draw("hist same");
    hist.data->Draw("p same");
    
    leg1->Clear();
    for(TH1F* sigvar : hist.sigvar) leg1->AddEntry( sigvar, sigvar->GetTitle(), "L" );
    leg1->AddEntry( hist.data, hist.data->GetTitle(), "P" );
    leg1->Draw();
    
    helper->DrawCMS();
    
    gPad->RedrawAxis();
    
    canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.data->GetName())+std::string("_sigvar.eps")).c_str(),"eps");
    canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.data->GetName())+std::string("_sigvar.png")).c_str(),"png");
    
  }
}
