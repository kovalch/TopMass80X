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
  std::string topBranchName(po::GetOption<std::string>("topBranchName"));
  int channelID_ = Helper::channelID();
  double lumi    = po::GetOption<double>("lumi");
  
  //
  // DEFINE HISTOGRAMS
  //

  // masses
  hists.push_back(MyHistogram("fitTop1Mass"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations", 75, 50, 400));
  hists.push_back(MyHistogram("fitTop1MassBest", "top.fitTop1[0].M()", "", ";m_{t}^{fit} [GeV]; Events"      , 75, 50, 400));
  hists.push_back(MyHistogram("recoTop1Mass"    , "top.recoTop1.M()"   , "", ";m_{t}^{reco} [GeV]; Permutations", 75, 50, 400));
  hists.push_back(MyHistogram("recoTop1MassBest", "top.recoTop1[0].M()", "", ";m_{t}^{reco} [GeV]; Events"      , 75, 50, 400));
  hists.push_back(MyHistogram("recoWAveMass"    , "(top.recoW1.M()+top.recoW2.M())/2.0"      , "", ";m_{W}^{reco} [GeV]; Permutations", 60, 65, 125));
  hists.push_back(MyHistogram("recoWAveMassBest", "(top.recoW1[0].M()+top.recoW2[0].M())/2.0", "", ";m_{W}^{reco} [GeV]; Events"      , 60, 65, 125));
  hists.push_back(MyHistogram("recoW1Mass"    , "top.recoW1.M()"   , "", ";m_{W}^{reco} [GeV]; Permutations", 60, 0, 300));
  hists.push_back(MyHistogram("recoW1MassBest", "top.recoW1[0].M()", "", ";m_{W}^{reco} [GeV]; Events"      , 60, 0, 300));

  hists.push_back(MyHistogram("fitTop1Mass_barrel", "top.fitTop1.M()", "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";m_{t}^{fit} [GeV], #eta^{j}<1.1; Permutations" , 75, 50, 400));
  hists.push_back(MyHistogram("recoW1Mass_barrel" , "top.recoW1.M()" , "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";m_{W}^{reco} [GeV], #eta^{j}<1.1; Permutations", 60,  0, 300));

  hists.push_back(MyHistogram("fitTop1Mass_vs_nVertex", "weight.nVertex", "top.fitTop1.M()", "", ";N_{Vertex}; m_{t}^{fit} [GeV]" , 10, 0, 50, 75, 50, 400));
  hists.push_back(MyHistogram("fitTop1Mass_vs_fitProb", "top.fitProb"   , "top.fitTop1.M()", "", ";P_{gof}; m_{t}^{fit} [GeV]"    , 10, 0,  1, 75, 50, 400));
  hists.push_back(MyHistogram("fitTop1Mass_vs_nJet", "jet.@jet.size()", "top.fitTop1.M()", "", ";N_{jet}; m_{t}^{fit} [GeV]", 3, 4, 7, 75, 50, 400));
  hists.push_back(MyHistogram("recoW1Mass_vs_nVertex" , "weight.nVertex", "top.recoW1.M()" , "", ";N_{Vertex}; m_{W}^{reco} [GeV]", 10, 0, 50, 58, 10, 300));
  hists.push_back(MyHistogram("recoW1Mass_vs_fitProb" , "top.fitProb"   , "top.recoW1.M()" , "", ";P_{gof}; m_{W}^{reco} [GeV]"   , 10, 0,  1, 58, 10, 300));
  hists.push_back(MyHistogram("recoW1Mass_vs_nJet", "jet.@jet.size()", "top.recoW1.M()", "", ";N_{jet}; m_{W}^{reco} [GeV]", 3, 4, 7, 58, 10, 300));

  // light pulls
  hists.push_back(MyHistogram("pullW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pull[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{q_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod2_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),-(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta())))))", "", ";#theta_{pull,ch}^{q_{2} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{q_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));

  hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_0",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi(),2)) < 1.",
     ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_1",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi(),2)) > 1.",
     ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_2",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi(),2)) > 2.",
    ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));

  // b pulls
  hists.push_back(MyHistogram("pullChargedB1_B2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))", "", ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{b_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB2_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitB2.Phi()),-(top.fitB1.Eta()-top.fitB2.Eta())))))", "", ";#theta_{pull,ch}^{b_{2} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{b_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));

  hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_0",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB2.Phi()-top.fitB1.Phi(),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) < 1.",
    ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_1",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB2.Phi()-top.fitB1.Phi(),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) > 1.",
    ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_2",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB2.Phi()-top.fitB1.Phi(),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) > 2.",
    ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));

  // mixed pulls
  hists.push_back(MyHistogram("pullChargedW1Prod1_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull}^{q_{1} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitB1.Phi()),-(top.fitW1Prod1.Eta()-top.fitB1.Eta())))))", "", ";#theta_{pull}^{b_{1} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));

  hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_0",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod1.Phi(),2)) < 1.",
     ";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_1",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod1.Phi(),2)) > 1.",
     ";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_2",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod1.Phi(),2)) > 2.",
     ";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));

  hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_0",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod2.Phi(),2)) < 1.",
     ";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_1",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod2.Phi(),2)) > 2.",
     ";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
  hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_2",
    "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
    "sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(top.fitB1.Phi()-top.fitW1Prod2.Phi(),2)) > 2.",
     ";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));

  // others
  hists.push_back(MyHistogram("fitProb"    , "top.fitProb"   , "", ";P_{gof}; Permutations", 50, 0, 1.0));
  hists.push_back(MyHistogram("fitProbBest", "top.fitProb[0]", "", ";P_{gof}; Events"      , 50, 0, 1.0));
  hists.push_back(MyHistogram("deltaRbb"    , "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))"            , "", ";#DeltaR_{b#bar{b}}; Permutations", 50, 1, 6));
  hists.push_back(MyHistogram("deltaRbbBest", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", "", ";#DeltaR_{b#bar{b}}; Events"      , 50, 1, 6));
  hists.push_back(MyHistogram("combinationType"    , "top.combinationType"   , "", ";Combination Type; Permutations", 20, -10, 10));
  hists.push_back(MyHistogram("combinationTypeBest", "top.combinationType[0]", "", ";Combination Type; Events"      , 20, -10, 10));

  // pts
  hists.push_back(MyHistogram("jet1Pt", "jet.jet[0].Pt()", "", ";p_{T}^{1} [GeV]; Events", 45, 0, 450));
  hists.push_back(MyHistogram("jet2Pt", "jet.jet[1].Pt()", "", ";p_{T}^{2} [GeV]; Events", 40, 0, 400));
  hists.push_back(MyHistogram("jet3Pt", "jet.jet[2].Pt()", "", ";p_{T}^{3} [GeV]; Events", 50, 0, 250));
  hists.push_back(MyHistogram("jet4Pt", "jet.jet[3].Pt()", "", ";p_{T}^{4} [GeV]; Events", 40, 0, 200));
  hists.push_back(MyHistogram("jet5Pt", "jet.jet[4].Pt()", "", ";p_{T}^{5} [GeV]; Events", 50, 0, 150));
  hists.push_back(MyHistogram("jet6Pt", "jet.jet[5].Pt()", "", ";p_{T}^{6} [GeV]; Events", 50, 0, 100));

  hists.push_back(MyHistogram("fitTop1Pt"    , "top.fitTop1.Pt()"   , "", ";p_{T,t}^{fit} [GeV]; Permutations", 40, 0, 400));
  hists.push_back(MyHistogram("fitTop1PtBest", "top.fitTop1[0].Pt()", "", ";p_{T,t}^{fit} [GeV]; Events"      , 40, 0, 400));
  hists.push_back(MyHistogram("fitTop1Pt_barrel", "top.fitTop1.Pt()", "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";p_{T,t}^{fit} [GeV], #eta^{j}<1.1; Permutations", 40, 0, 400));
  hists.push_back(MyHistogram("fitRlb", "(top.fitB1.Pt()+top.fitB2.Pt())/(top.fitW1Prod1.Pt()+top.fitW1Prod2.Pt())", "", ";R_{lb}^{fit}; Permutations", 40, 0, 4));
  hists.push_back(MyHistogram("recoRlb", "(top.recoB1.Pt()+top.recoB2.Pt())/(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())", "", ";R_{lb}^{reco}; Permutations", 40, 0, 4));

  hists.push_back(MyHistogram("fitB1Pt", "top.fitB1.Pt()", "", ";p_{T}^{B1} [GeV]; Permutations", 40, 0, 200));
  hists.push_back(MyHistogram("fitB2Pt", "top.fitB2.Pt()", "", ";p_{T}^{B2} [GeV]; Permutations", 40, 0, 200));
  hists.push_back(MyHistogram("fitW1Prod1Pt", "top.fitW1Prod1.Pt()", "", ";p_{T}^{W1,1} [GeV]; Permutations", 40, 0, 200));
  hists.push_back(MyHistogram("fitW1Prod2Pt", "top.fitW1Prod2.Pt()", "", ";p_{T}^{W1,2} [GeV]; Permutations", 40, 0, 200));
  hists.push_back(MyHistogram("fitW2Prod1Pt", "top.fitW2Prod1.Pt()", "", ";p_{T}^{W2,1} [GeV]; Permutations", 40, 0, 200));
  hists.push_back(MyHistogram("fitW2Prod2Pt", "top.fitW2Prod2.Pt()", "", ";p_{T}^{W2,2} [GeV]; Permutations", 40, 0, 200));

  // lepton+jets
  hists.push_back(MyHistogram("leptonPt", "top.fitW2Prod1.Pt()", "", ";p_{T}^{lepton} [GeV]; Events", 40, 0, 200));
  hists.push_back(MyHistogram("leptonEta", "top.fitW2Prod1.Eta()", "", ";#eta^{lepton}; Events", 25, -2.5, 2.5));

  // jet details
  hists.push_back(MyHistogram("nChargedHadronsBJetT", "jet.nChargedHadrons", "jet.bTagCSV>0.898", ";N_{ch}; B Jets", 40, 0, 200));
  hists.push_back(MyHistogram("nChargedHadronsBJetM", "jet.nChargedHadrons", "jet.bTagCSV>0.679", ";N_{ch}; B Jets", 40, 0, 200));
  hists.push_back(MyHistogram("nChargedHadronsLJetT", "jet.nChargedHadrons", "jet.bTagCSV<0.898", ";N_{ch}; Non B Jets", 40, 0, 200));
  hists.push_back(MyHistogram("nChargedHadronsLJetM", "jet.nChargedHadrons", "jet.bTagCSV<0.679", ";N_{ch}; Non B Jets", 40, 0, 200));

  // event observables
  hists.push_back(MyHistogram("nJet", "jet.@jet.size()", "", ";N_{jet}; Events", 15, 0, 15));
  hists.push_back(MyHistogram("nVertex", "weight.nVertex", "", ";N_{vertex}; Events", 50, 0, 50));

  MyBRegVarInfo helperMyBRegVarInfo;
  for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
	hists.push_back(MyHistogram("B1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
    hists.push_back(MyHistogram("B1HighGBR"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", "BRegJet.BRegGBRTrainResult["+topBranchName+"recoJetIdxB1]>1.6",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
    hists.push_back(MyHistogram("W1Prod1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxW1Prod1]", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
    //2D
    hists.push_back(MyHistogram("B1TopMassVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", topBranchName+"fitTop1.M()", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; m_{t}^{fit} [GeV]",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 50, 400 ));
    hists.push_back(MyHistogram("B1GBRTrainFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", "BRegJet.BRegGBRTrainResult["+topBranchName+"recoJetIdxB1]", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; GBRFactor",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0., 2. ));
    hists.push_back(MyHistogram("B1recoRlbVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", "(top.recoB1.Pt()+top.recoB2.Pt())/(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";R_{lb}^{reco}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 4));
    hists.push_back(MyHistogram("B1TrueCorrFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"["+topBranchName+"recoJetIdxB1]", "BRegJet.genJetPt["+topBranchName+"recoJetIdxB1]/BRegJet.jetPtCorr["+topBranchName+"recoJetIdxB1]", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{gen}/p_{T}^{corr}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));

  }


  //
  // DEFINE DATASETS
  //
  
  // Alljets channel
  if (channelID_ == Helper::kAllJets) {
    samples.push_back(MySample("Data", "MJP12*v1_data", kData, kBlack));
    samples.push_back(MySample("t#bar{t}", "Z2_S12_ABS_JES_100_172_5_sig", kSig, kRed+1));
    samples.push_back(MySample("QCD", "QCDMixing_MJPS12*v1_data", kBkg, kYellow, 1, 3.12));

    /* Signal variations (UE Tune)
    samples.push_back(MySample("t#bar{t}, Z2*"     , "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kRed+1    , 1));
    samples.push_back(MySample("t#bar{t}, P11"     , "Z2_S12_P11_sig"              , kSigVar, kMagenta+1, 9));
    samples.push_back(MySample("t#bar{t}, P11mpiHi", "Z2_S12_P11mpiHi_sig"         , kSigVar, kBlue+1   , 3));
    samples.push_back(MySample("t#bar{t}, P11TeV"  , "Z2_S12_P11TeV_sig"           , kSigVar, kGreen+4  , 4));
    samples.push_back(MySample("t#bar{t}, P11noCR" , "Z2_S12_P11NoCR_sig"          , kSigVar, kCyan+1   , 2));
    //*/
    /* Signal variations (MC Generator)
    samples.push_back(MySample("t#bar{t}, Z2*"          , "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kRed+1  , 1));
    samples.push_back(MySample("t#bar{t}, Powheg+Pythia", "Z2_S12_POWHEG_sig"           , kSigVar, kGreen+1, 9));
    samples.push_back(MySample("t#bar{t}, Powheg+Herwig", "Z2_S12_POWHER_sig"           , kSigVar, kBlue+1 , 7));
    samples.push_back(MySample("t#bar{t}, MC@NLO"       , "Z2_S12_MCNLO_sig"            , kSigVar, kCyan+1 , 2));
    /*/
    /* Signal variations (JES)
    samples.push_back(MySample("t#bar{t}, JES = 0.96", "Z2_S12_ABS_JES_096_172_5_sig", kSigVar, kRed+1    , 1));
    samples.push_back(MySample("t#bar{t}, JES = 0.98", "Z2_S12_ABS_JES_098_172_5_sig", kSigVar, kMagenta+1, 1));
    samples.push_back(MySample("t#bar{t}, JES = 1.00", "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kBlue+1   , 1));
    samples.push_back(MySample("t#bar{t}, JES = 1.02", "Z2_S12_ABS_JES_102_172_5_sig", kSigVar, kCyan+1   , 1));
    samples.push_back(MySample("t#bar{t}, JES = 1.04", "Z2_S12_ABS_JES_104_172_5_sig", kSigVar, kGreen+1  , 1));
    //*/
    /* Signal variations (MC Modelling)
    samples.push_back(MySample("t#bar{t}, , "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kRed+1    , 1));
    samples.push_back(MySample("t#bar{t}, , "Z2_S12_Scale_Up_sig"         , kSigVar, kMagenta+1, 9));
    samples.push_back(MySample("t#bar{t}, , "Z2_S12_Scale_Down_sig"       , kSigVar, kBlue+1   , 7));
    samples.push_back(MySample("t#bar{t}, , "Z2_S12_Matching_Up_sig"      , kSigVar, kCyan+1   , 2));
    samples.push_back(MySample("t#bar{t}, , "Z2_S12_Matching_Down_sig"    , kSigVar, kGreen+1  , 3));
    //*/
  }
  
  // Lepton+jets channel
  else {
    samples.push_back(MySample("Data", "Run2012", kData, kBlack));
    //    samples.push_back(MySample("t#bar{t} w.o. b-regression", "Summer12_TTJets1725_1.00", kData, kBlack));
    samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kSig, kRed+1, 1, lumi/1000.));
    samples.push_back(MySample("Z+Jets", "Summer12_ZJets", kBkg, kAzure-2, 1, lumi/1000.));
    samples.push_back(MySample("W+Jets", "Summer12_WJets", kBkg, kGreen-3, 1, lumi/1000.));
    samples.push_back(MySample("single top", "Summer12_singleTop", kBkg, kMagenta, 1, lumi/1000.));
    
    samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi/1000.));
    //* Signal variations
    samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_MGDecays", kSigVar, kRed+1, 1, lumi/1000.*9./4.));
    samples.push_back(MySample("t#bar{t}, P11", "Summer12_TTJets1725_MGDecays_P11", kSigVar, kMagenta+1, 9, lumi/1000.*9./4.));
    samples.push_back(MySample("t#bar{t}, P11noCR", "Summer12_TTJets1725_MGDecays_P11noCR", kSigVar, kCyan+1, 2, lumi/1000.*9./4.));
    //*/
    /* Signal variations (Powheg)
    samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, Powheg+Pythia", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, Powheg+Herwig", "Summer12_TTJets1725_powheg_herwig", kSigVar, kBlue+1, 7, lumi/1000.));
    //*/
    /* Signal variations (JES)
    samples.push_back(MySample("t#bar{t}, JES = 0.96", "Summer12_TTJets1725_0.96", kSigVar, kRed+1, 1, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, JES = 0.98", "Summer12_TTJets1725_0.98", kSigVar, kMagenta+1, 1, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, JES = 1.00", "Summer12_TTJets1725_1.00", kSigVar, kBlue+1, 1, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, JES = 1.02", "Summer12_TTJets1725_1.02", kSigVar, kCyan+1, 1, lumi/1000.));
    samples.push_back(MySample("t#bar{t}, JES = 1.04", "Summer12_TTJets1725_1.04", kSigVar, kGreen+1, 1, lumi/1000.));
    //*/
  }
  
  //
  // FILL HISTOGRAMS
  //

  {
    int bkgCounter = -1;
    int sigVarCounter = -1;
    // Loop over all samples
    for(MySample& sample : samples){
      if(sample.type == kBkg    ) ++bkgCounter;
      if(sample.type == kSigVar ) ++sigVarCounter;
      // Get sample
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
      std::cout << "Adding " << nFiles << " files for " << sample.name << " (" << sample.file << "), type: " << sample.type << std::endl;

      // Soft-initialize histograms and add sample
      for(MyHistogram& hist : hists) {
        if(hist.Dimension() == -1) continue;
        hist.Init(chain, topBranchName);
        if     (sample.type == kSig    ) hist.AddSignal(&sample);
        else if(sample.type == kBkg    ) hist.AddBackground(&sample);
        else if(sample.type == kSigVar ) hist.AddSignalVariation(&sample);
      }

      // Initialize weight and selection formulas
      TTreeFormula weight("weight",  po::GetOption<std::string>("weight").c_str(), chain);
      TTreeFormula sel   ("sel"   ,  po::GetOption<std::string>("analysisConfig.selection").c_str(), chain);
      TTreeFormula selCP ("selCP" , (po::GetOption<std::string>("analysisConfig.selection")
				 +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionCP")).c_str(), chain);
      TTreeFormula selWP ("selWP" , (po::GetOption<std::string>("analysisConfig.selection")
				 +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionWP")).c_str(), chain);
      TTreeFormula selUN ("selUN" , (po::GetOption<std::string>("analysisConfig.selection")
                 +std::string(" & ")+po::GetOption<std::string>("analysisConfig.selectionUN")).c_str(), chain);

      //  Loop over all events
      for(int i = 0; ; ++i){
        long entry = chain->LoadTree(i);
        if(entry  < 0) break;
        if(entry == 0){
          for(auto& hist : hists) {
            hist.varx->UpdateFormulaLeaves();
            hist.vary->UpdateFormulaLeaves();
            if (hist.selection.size() > 0) hist.sel->UpdateFormulaLeaves();
          }
          weight.UpdateFormulaLeaves();
          sel   .UpdateFormulaLeaves();
          selCP .UpdateFormulaLeaves();
          selWP .UpdateFormulaLeaves();
          selUN .UpdateFormulaLeaves();
        }
        // Skip event if base selection is not met
        if(!weight.GetNdata()) continue;
        if(!sel   .GetNdata()) continue;
        selCP.GetNdata(); selWP.GetNdata(); selUN.GetNdata();

        // Loop over permutations
        for(int j = 0, l = sel.GetNdata(); j < l; ++j){
          if(!sel.EvalInstance(j)) continue;
          // Loop over samples
          for(MyHistogram& hist : hists){
        	if(!hist.varx->GetNdata()) continue;
            if(hist.Dimension() == -1) continue;
            if(hist.varx->GetNdata()<=j) continue;
            if(hist.Dimension() == 1){
              if(!hist.vary->GetNdata()) continue;
            }
            else if(hist.Dimension() == 1){
              if(hist.vary->GetNdata()<=j) continue;
            }
	        else if(hist.Dimension() == 2){
	       	  if(!hist.vary->GetNdata()) continue;
	       	}
            // Skip permutation if individual selection is not met
            if (hist.selection.size() > 0) {
              if(!hist.sel->GetNdata()) continue;
              if(!hist.sel->EvalInstance(j)) continue;
            }
            // Fill according to sample and histogram type
            // Fill data
            if     (sample.type == kData && hist.Dimension() == 1) {
              hist.Data1D()->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kData && hist.Dimension() == 2) {
              hist.Data2D()->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            // Fill signal
            else if(sample.type == kSig && hist.Dimension() == 1) {
              if     (selCP.EvalInstance(j)) hist.Sig1D().at(2)->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selWP.EvalInstance(j)) hist.Sig1D().at(1)->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selUN.EvalInstance(j)) hist.Sig1D().at(0)->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kSig && hist.Dimension() == 2) {
              if     (selCP.EvalInstance(j)) hist.Sig2D().at(2)->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selWP.EvalInstance(j)) hist.Sig2D().at(1)->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              else if(selUN.EvalInstance(j)) hist.Sig2D().at(0)->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            // Fill background
            else if(sample.type == kBkg  && hist.Dimension() == 1) {
              hist.Bkg1D()[bkgCounter]->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              // Add background also to histogram for signal variation
              // TODO: Background normalization in signal variations not correct for AllJets
              for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kBkg  && hist.Dimension() == 2) {
              hist.Bkg2D()[bkgCounter]->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
              // Add background also to histogram for signal variation
              // TODO: Background normalization in signal variations not correct for AllJets
              for(TH2F* sigvar : hist.Sigvar2D()) sigvar->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            }
            // Fill signal variation
            else if(sample.type == kSigVar && hist.Dimension() == 1) hist.Sigvar1D()[sigVarCounter]->Fill(hist.varx->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
            else if(sample.type == kSigVar && hist.Dimension() == 2) hist.Sigvar2D()[sigVarCounter]->Fill(hist.varx->EvalInstance(j), hist.vary->EvalInstance(j), weight.EvalInstance(j)*sample.scale);
          }
        }
      }
      for(MyHistogram& hist : hists) {
        hist.varx->Clear();
        hist.vary->Clear();
      }
    }
  }
  
  //
  // DRAW CONTROL PLOTS
  //
  
  bool firstHist = true;
  for(MyHistogram& hist : hists){
    if(hist.Dimension() == -1) continue;

    // Show event yields for first histogram
    int bins = hist.Data1D()->GetNbinsX()+1;
    if (firstHist) std::cout << "Yields" << std::endl;
    double integralD = hist.Data1D()->Integral(0,bins);
    if (firstHist) std::cout << "Data:       " << integralD << std::endl;
    double integralS = 0;
    for(TH1F* sig : hist.Sig1D()) {
      if (firstHist) std::cout << "  " << sig->GetTitle() << ": " << sig->Integral(0,bins) << std::endl;
      integralS += sig->Integral(0,bins);
    }
    if (firstHist) std::cout << "Signal:     " << integralS << std::endl;
    double integralB = 0;
    for(TH1F* bkg : hist.Bkg1D()) {
      if (firstHist) std::cout << "  " << bkg->GetTitle() << ": " << bkg->Integral(0,bins) << std::endl;
      integralB += bkg->Integral(0,bins);
    }
    if (firstHist) std::cout << "Background: " << integralB << std::endl;
    if (firstHist) std::cout << "MC Total:   " << integralS+integralB << std::endl;
    firstHist = false;
    
    // Alljets: Normalize to data, use fixed signal fraction
    if (channelID_ == Helper::kAllJets) {
      double fSig = po::GetOption<double>("templates.fSig");
      for(TH1F* sig    : hist.Sig1D())    sig->Scale(    fSig *integralD/integralS);
      for(TH1F* bkg    : hist.Bkg1D())    bkg->Scale((1.-fSig)*integralD/integralB);
      for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Scale(integralD/sigvar->Integral(0,bins));
    }
    // Lepton+jets: Normalize to data
    else {
      for(TH1F* sig    : hist.Sig1D())    sig->Scale(integralD/(integralS+integralB));
      for(TH1F* bkg    : hist.Bkg1D())    bkg->Scale(integralD/(integralS+integralB));
      for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Scale(integralD/sigvar->Integral(0,bins));
      
      // Double bin error due to correlations
      // (Many observables do not change for different neutrino solutions or b assignments)
      for (int i = 0; i < hist.Data1D()->GetNbinsX(); ++i) {
        hist.Data1D()->SetBinError(i, hist.Data1D()->GetBinError(i) * 2);
      }
    }
    
    // Draw 2D plots
    if(hist.Dimension() == 2){
      std::cout << "doing 2D plots" << std::endl;

      TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
      canv->cd();
      canv->SetRightMargin(0.15);


      std::vector<TH2F*> collectAll2D;
      std::vector<TH1*> collectAll2DProfiles;
      collectAll2D.push_back(hist.Data2D());
      for(TH2F* sig : hist.Sig2D()) collectAll2D.push_back(sig);
      for(TH2F* bkg : hist.Bkg2D()) collectAll2D.push_back(bkg);

      for (size_t h_i=0;h_i<collectAll2D.size();++h_i){
    	  collectAll2D.at(h_i)->SetMarkerStyle(20+h_i);
    	  collectAll2D.at(h_i)->SetLineColor(collectAll2D.at(h_i)->GetFillColor());
    	  collectAll2D.at(h_i)->SetMarkerColor(collectAll2D.at(h_i)->GetFillColor());
    	  if(h_i==0){
    	    collectAll2D.at(h_i)->SetMarkerColor(1);
    		collectAll2D.at(h_i)->SetLineColor(1);
    	  }
    	  collectAll2DProfiles.push_back(collectAll2D.at(h_i)->ProfileX(collectAll2D.at(h_i)->GetName()+(TString)"_pfx"));
    	  collectAll2DProfiles.at(h_i)->GetYaxis()->SetTitle(collectAll2D.at(h_i)->GetYaxis()->GetTitle());
    	  std::cout << "collectAll2D.at(h_i)->GetName()+(TString)\"_pfx\" " << collectAll2D.at(h_i)->GetName()+(TString)"_pfx" << std::endl;
    	  collectAll2D.at(h_i)->Draw("colz");
    	  collectAll2DProfiles.at(h_i)->Draw("same");

    	  TLegend* legInd = new TLegend(0.25, 0.85, 0.55, 0.925);
    	  legInd->SetTextSize(0.03);
    	  legInd->SetFillStyle(0);
    	  legInd->SetBorderSize(0);
    	  legInd->AddEntry( collectAll2D.at(h_i), collectAll2D.at(h_i)->GetTitle(), "LP" );
    	  legInd->Draw();

    	  helper->DrawCMS();
    	  gPad->RedrawAxis();
    	  std::cout << HelperFunctions::cleanedName(collectAll2D.at(h_i)->GetTitle()) << std::endl;
    	  canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data2D()->GetName())+std::string("_Ind2D_")+HelperFunctions::cleanedName(collectAll2D.at(h_i)->GetTitle())+std::string(".eps")).c_str(),"eps");
      }


      HelperFunctions::setCommonYRange(collectAll2DProfiles,0.35);
      
      canv->SetRightMargin(0.05);

      for (size_t h_i=0;h_i<collectAll2DProfiles.size();++h_i){
    	  if(h_i==0)collectAll2DProfiles.at(h_i)->Draw();
    	  else collectAll2DProfiles.at(h_i)->Draw("SAME");
      }

      TLegend* leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
      leg0->SetTextSize(0.03);
      leg0->SetFillStyle(0);
      leg0->SetBorderSize(0);
      for(TH2F* sig : hist.Sig2D()){
        leg0->AddEntry( sig, sig->GetTitle(), "LP" );
      }
      leg0->Draw();

      TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
      leg1->SetTextSize(0.03);
      leg1->SetFillStyle(0);
      leg1->SetBorderSize(0);
      for(TH2F* bkg : hist.Bkg2D()) leg1->AddEntry( bkg, bkg->GetTitle(), "LP" );
      leg1->AddEntry( hist.Data2D(), hist.Data2D()->GetTitle(), "LP" );
      leg1->Draw();

      helper->DrawCMS();
      gPad->RedrawAxis();


      canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data2D()->GetName())+std::string(".eps")).c_str(),"eps");
      canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data2D()->GetName())+std::string(".png")).c_str(),"png");


      collectAll2D.clear();
      for(TH1* prof : collectAll2DProfiles) delete prof;//delete profiles
      collectAll2DProfiles.clear();

      // Draw signal variation plots only if there are signal variations
      if(hist.Sigvar2D().size()){
        collectAll2D.push_back(hist.Data2D());
        for(TH2F* sigvar : hist.Sigvar2D()) collectAll2D.push_back(sigvar);


        for (size_t h_i=0;h_i<collectAll2D.size();++h_i){
          std::cout << " title " << collectAll2D.at(h_i)->GetTitle() << " " << collectAll2D.at(h_i)->GetFillColor() << std::endl;
          collectAll2D.at(h_i)->SetMarkerStyle(20+h_i);
          //	collectAll2D.at(h_i)->SetLineColor(collectAll2D.at(h_i)->GetFillColor());
          collectAll2D.at(h_i)->SetMarkerColor(collectAll2D.at(h_i)->GetLineColor());
          if(h_i==0){
            collectAll2D.at(h_i)->SetMarkerColor(1);
            collectAll2D.at(h_i)->SetLineColor(1);
          }
          collectAll2DProfiles.push_back(collectAll2D.at(h_i)->ProfileX(collectAll2D.at(h_i)->GetName()+(TString)"_pfx"));
          collectAll2DProfiles.at(h_i)->GetYaxis()->SetTitle(collectAll2D.at(h_i)->GetYaxis()->GetTitle());
        }

        HelperFunctions::setCommonYRange(collectAll2DProfiles,0.35);

        for (size_t h_i=0;h_i<collectAll2DProfiles.size();++h_i){
          if(h_i==0){
            collectAll2DProfiles.at(h_i)->Draw();
          }
          else {
            collectAll2DProfiles.at(h_i)->Draw("SAME");
          }
        }

        leg1->Clear();
        for(TH2F* sigvar : hist.Sigvar2D()) leg1->AddEntry( sigvar, sigvar->GetTitle(), "LP" );
        leg1->AddEntry( hist.Data2D(), hist.Data2D()->GetTitle(), "P" );
        leg1->Draw();

        helper->DrawCMS();

        gPad->RedrawAxis();

        canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.eps")).c_str(),"eps");
        canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.png")).c_str(),"png");
      }//end if 2D signal variation plots
    }
    else if(hist.Dimension() == 1){

      // Draw control plot
      TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
      canv->cd();

      hist.Data1D()->GetYaxis()->SetRangeUser(1, hist.Data1D()->GetMaximum()*1.5);
      hist.Data1D()->Draw("p");
      THStack* stack = new THStack("stack", "");
      //    std::cout << "building stack hist.bkg.size() " << hist.bkg.size() << " hist.Bkg1D().size() " <<  hist.Bkg1D().size() << std::endl;
      //    std::cout << hist.Bkg1D().at(0) << " " << hist.bkg.at(0) << std::endl;
      //    std::cout << hist.Bkg1D().at(1) << " " << hist.bkg.at(1) << std::endl;
      //    std::cout << hist.Bkg1D().at(2) << " " << hist.bkg.at(2) << std::endl;


      //unfortunately PYTHON-style for loop breaks when doing boost::adaptors::reverse(hist.Bkg1D()) directly (pointers get mixed up)
      std::vector <TH1F*> tempHistBkg1D = hist.Bkg1D();
      for(TH1F* bkg : boost::adaptors::reverse(tempHistBkg1D)) stack->Add(bkg);
      std::vector <TH1F*> tempHistSig1D = hist.Sig1D();
      for(TH1* sig : boost::adaptors::reverse(tempHistSig1D)) stack->Add(sig);

      stack     ->Draw("hist same");
      hist.Data1D()->Draw("p same");

      TLegend* leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
      leg0->SetTextSize(0.03);
      leg0->SetFillStyle(0);
      leg0->SetBorderSize(0);
      for(TH1F* sig : hist.Sig1D()) leg0->AddEntry( sig, sig->GetTitle(), "F" );
      leg0->Draw();

      TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
      leg1->SetTextSize(0.03);
      leg1->SetFillStyle(0);
      leg1->SetBorderSize(0);
      for(TH1F* bkg : hist.Bkg1D()) leg1->AddEntry( bkg, bkg->GetTitle(), "F" );
      leg1->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "P" );
      leg1->Draw();

      helper->DrawCMS();

      gPad->RedrawAxis();

      canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string(".eps")).c_str(),"eps");
      canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string(".png")).c_str(),"png");


      //draw together with ratio underneath
      TH1* ratioToTHStack = HelperFunctions::createRatioPlot((TH1 *)hist.Data1D(), ((TH1 *)(stack->GetStack()->Last())), "Data/MC");
      TCanvas* canvWRatio = new TCanvas("cControlPlotsWRatio", "cControlPlotsWRatio", 600, 600);
      canvWRatio->Range(0,0,1,1);
      canvWRatio->cd();

      TPad *pad2 = new TPad("pad2","pad2",0,0,1,1);
      pad2->SetTopMargin(0.71);
      pad2->Draw();
      pad2->cd();
      ratioToTHStack->Draw();
      ratioToTHStack->GetYaxis()->SetRangeUser(0.49,1.51);
      ratioToTHStack->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
      ratioToTHStack->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X")*0.7);
      ratioToTHStack->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X")*0.7);
      ratioToTHStack->GetYaxis()->SetNdivisions(205);
      ratioToTHStack->GetYaxis()->CenterTitle();
      ratioToTHStack->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y")*0.7);
      ratioToTHStack->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y")*0.7);
      ratioToTHStack->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()/0.7);


      canvWRatio->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
      pad1->SetBottomMargin(0.0);
      pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      pad1->Draw();
      pad1->cd();
      hist.Data1D()->Draw("p");
      hist.Data1D()->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
      hist.Data1D()->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y"));
      stack     ->Draw("hist same");
      hist.Data1D()->Draw("p same");
      leg0->Draw();
      leg1->Draw();

      canvWRatio->cd();

      helper->DrawCMS();
      gPad->RedrawAxis();

      canvWRatio->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_Ratio.eps")).c_str(),"eps");
      canvWRatio->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_Ratio.png")).c_str(),"png");

      // Draw signal variation plot only if there are signal variations
      if(hist.Sigvar1D().size()){
        canv->cd();
        canv->Clear();

        hist.Data1D()->GetYaxis()->UnZoom();
        hist.Data1D()->GetYaxis()->SetRangeUser(hist.Data1D()->GetMinimum()*0.9, hist.Data1D()->GetMaximum()*1.15);
        hist.Data1D()->Draw("p");
        for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Draw("hist same");
        hist.Data1D()->Draw("p same");

        leg1->Clear();
        for(TH1F* sigvar : hist.Sigvar1D()) leg1->AddEntry( sigvar, sigvar->GetTitle(), "L" );
        leg1->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "P" );
        leg1->Draw();

        helper->DrawCMS();

        gPad->RedrawAxis();

        canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.eps")).c_str(),"eps");
        canv->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.png")).c_str(),"png");


        //ratio plots
        canvWRatio->cd();
        //      canvWRatio->Clear();
        std::vector<TH1*> collectRatios;
        for(TH1F* sigvar : hist.Sigvar1D()){
          collectRatios.push_back(HelperFunctions::createRatioPlot((TH1 *)hist.Data1D(), ((TH1 *)(sigvar)), "Data/MC"));
          collectRatios.back()->SetMarkerColor(sigvar->GetLineColor());
          collectRatios.back()->SetLineColor(sigvar->GetLineColor());
          collectRatios.back()->SetLineStyle(sigvar->GetLineStyle());
        }

        pad2->cd();
        pad2->Draw();
        for (size_t h_i=0;h_i<collectRatios.size();++h_i){
          if(h_i==0){
            collectRatios.at(h_i)->Draw("hist");
            collectRatios.at(h_i)->GetYaxis()->SetRangeUser(0.49,1.51);
            collectRatios.at(h_i)->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
            collectRatios.at(h_i)->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X")*0.7);
            collectRatios.at(h_i)->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X")*0.7);
            collectRatios.at(h_i)->GetYaxis()->SetNdivisions(205);
            collectRatios.at(h_i)->GetYaxis()->CenterTitle();
            collectRatios.at(h_i)->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y")*0.7);
            collectRatios.at(h_i)->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y")*0.7);
            collectRatios.at(h_i)->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()/0.7);
          }
          else collectRatios.at(h_i)->Draw("hist same");
        }
        pad1->Draw();
        pad1->cd();
        hist.Data1D()->GetYaxis()->UnZoom();
        hist.Data1D()->GetYaxis()->SetRangeUser(hist.Data1D()->GetMinimum()*0.9, hist.Data1D()->GetMaximum()*1.15);
        hist.Data1D()->Draw("p");
        for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Draw("hist same");
        hist.Data1D()->Draw("p same");
        leg1->Draw();

        canvWRatio->cd();
        helper->DrawCMS();
        gPad->RedrawAxis();


        canvWRatio->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar_Ratio.eps")).c_str(),"eps");
        canvWRatio->Print((std::string("plot/controlplots/")+channel+std::string("/")+channel+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar_Ratio.png")).c_str(),"png");
      }//end if 1D signal variation plots
    }//end 1D plots
    else{
      std::cout << "The Histogram *" << hist.Data1D()->GetName() << "* has an invalid dimension *" << hist.Dimension() << "*! Supported are only 1 and 2!" << std::endl;
    }
  }
}
