#include "controlPlots.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TDirectory.h"
//#include "TLegend.h"
//#include "TROOT.h"
//#include "TString.h"
//#include "TGaxis.h"
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
#include "TLegendEntry.h"
#include "TFile.h"

#include "ProgramOptionsReader.h"
#include "Helper.h"
#include "CMS_lumi.h"

#include "../../TopEventTree/interface/MyMa.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/range/adaptor/reversed.hpp>


//FIXME doesnt create its "plot/controlplots" directory, add that

typedef ProgramOptionsReader po;

TopMassControlPlots::TopMassControlPlots() :
      path_         (po::GetOption<std::string>("analysisConfig.samplePath")),
      outPath_      (po::GetOption<std::string>("outPath")),
      channel_      (po::GetOption<std::string>("channel")),
      topBranchName_(po::GetOption<std::string>("topBranchName")),
      lumi_         (po::GetOption<double>("lumi"))
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
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(1.5);

  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  int channelID = Helper::channelID();
  
  std::map<std::string,bool> plotSelectedForPlotting;
  std::vector<std::string> selectPlotsToDraw  = helper->readParametersString("analysisConfig.plotsToDraw");
  std::cout << "These plots will be drawn: " << std::endl;
  std::cout << "===========================" << std::endl;
  for(std::string plotToDraw : selectPlotsToDraw){
	  std::cout << plotToDraw << std::endl;
	  plotSelectedForPlotting[plotToDraw]=true;
  }
  std::cout << "===========================" << std::endl;
  std::cout << "Extra label for mass plots: " << po::GetOption<std::string>("analysisConfig.extraLabel") << std::endl;
  std::cout << "===========================" << std::endl;
  

  //
  // DEFINE HISTOGRAMS
  //
  
  hists.push_back(MyHistogram("leptonFlavour", "top.leptonFlavour[0]"   , "", ";Lepton flavour; Events", 40, -20, 20));
  hists.push_back(MyHistogram("combinationType"    , "top.combinationType"   , "", ";Combination Type; Permutations", 20, -10, 10));

  // others (do these first to get correct event yield)
  if(plotSelectedForPlotting.find("test")!=plotSelectedForPlotting.end()){
    //hists.push_back(MyHistogram("leptonFlavour", "top.leptonFlavour[0]"   , "", ";Lepton flavour; Events", 40, -20, 20));
    hists.push_back(MyHistogram("Nbjet", "Sum$(jet.jet.Pt()>30 & jet.bTagCSV>0.679)"   , "", ";b-jet multiplicity; Events", 6, -0.5, 5.5));
    hists.push_back(MyHistogram("fitTop1Mass"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
  }
  
  if(plotSelectedForPlotting.find("ExtraPlotsFitCombTypeEtc")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitProb"    , "top.fitProb"   , "", ";P_{gof}; Permutations / 0.02", 50, 0, 1.0));
    hists.push_back(MyHistogram("fitProbLogY", "top.fitProb"   , "", ";P_{gof}; Permutations / 0.02", 50, 0, 1.0));
    hists.back().ConfigureExtraOptions(false, "1.", false, "", false, true);
    hists.push_back(MyHistogram("fitProbNorm", "top.fitProb"   , "", ";P_{gof}; Permutations", 50, 0, 1.0));
    hists.back().ConfigureMoreExtraOptions(true);
    hists.push_back(MyHistogram("fitProbLow" , "top.fitProb"   , "", ";P_{gof}; Permutations / 0.2",  1, 0, 0.2));
    hists.push_back(MyHistogram("fitProbBest", "top.fitProb[0]", "", ";P_{gof}; Events / 0.02"      , 50, 0, 1.0));
    hists.push_back(MyHistogram("fitProbBestLogY", "top.fitProb[0]", "", ";P_{gof}; Events / 0.02"  , 50, 0, 1.0));
    hists.back().ConfigureExtraOptions(false, "1.", false, "", false, true);
    hists.push_back(MyHistogram("fitChi2"    , "top.fitChi2"   , "", ";#chi^{2}; Permutations / 1", 50, 0, 50));
    hists.push_back(MyHistogram("deltaRbb"    , "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))"            , "", ";#DeltaR_{b#bar{b}}; Permutations / 0.1", 50, 1, 6));
    hists.push_back(MyHistogram("deltaRbbBest", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", "", ";#DeltaR_{b#bar{b}}; Events / 0.1"      , 50, 1, 6));
    hists.push_back(MyHistogram("deltaRbbLow"    , "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))"            , "", ";#DeltaR_{b#bar{b}}; Permutations / 0.1", 60, 0, 6));
    hists.push_back(MyHistogram("deltaRbbBestLow", "sqrt(pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2))", "", ";#DeltaR_{b#bar{b}}; Events / 0.1"      , 60, 0, 6));
    hists.push_back(MyHistogram("deltaRqq"    , "sqrt(pow(top.recoW1Prod1.Eta()-top.recoW1Prod2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.recoW1Prod1.Phi()-top.recoW1Prod2.Phi()),2))", "", ";#DeltaR_{q#bar{q}}; Permutations / 0.1", 40, 0, 4));
    hists.push_back(MyHistogram("combinationType"    , "top.combinationType"   , "", ";Combination Type; Permutations", 20, -10, 10));
    hists.push_back(MyHistogram("combinationTypeBest", "top.combinationType[0]", "", ";Combination Type; Events"      , 20, -10, 10));
  }
  
  // masses
  if(plotSelectedForPlotting.find("BasicMasses")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTop1Mass"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
    hists.push_back(MyHistogram("fitTop1MassBest", "top.fitTop1[0].M()", "", ";m_{t}^{fit} [GeV]; Events / 5 GeV"      , 70, 50, 400));
    //hists.back().SetFitGaussToCore();
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
    hists.push_back(MyHistogram("fitTop1MassPeak"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations / 2 GeV", 45, 125, 215));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassPeakBest", "top.fitTop1[0].M()", "", ";m_{t}^{fit} [GeV]; Events / 2 GeV"      , 45, 125, 215));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("recoTop1Mass"    , "top.recoTop1.M()"   , "", ";m_{t}^{reco} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
    hists.push_back(MyHistogram("recoTop1MassBest", "top.recoTop1[0].M()", "", ";m_{t}^{reco} [GeV]; Events / 5 GeV"      , 70, 50, 400));
    hists.push_back(MyHistogram("recoWAveMass"    , "(top.recoW1.M()+top.recoW2.M())/2.0"      , "", ";m_{W}^{reco} [GeV]; Permutations / 1 GeV", 60, 65, 125));
    hists.push_back(MyHistogram("recoWAveMassBest", "(top.recoW1[0].M()+top.recoW2[0].M())/2.0", "", ";m_{W}^{reco} [GeV]; Events / 1 GeV"      , 60, 65, 125));
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
    hists.push_back(MyHistogram("recoW1Mass"    , "top.recoW1.M()"   , "", ";m_{W}^{reco} [GeV]; Permutations / 5 GeV", 60, 0, 300));
    hists.back().SetExtraLabel(po::GetOption<std::string>("analysisConfig.extraLabel"));
    hists.push_back(MyHistogram("recoW1MassBest", "top.recoW1[0].M()", "", ";m_{W}^{reco} [GeV]; Events / 5 GeV"      , 60, 0, 300));
  }
  if(plotSelectedForPlotting.find("ExtraMasses")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTop1Mass_barrel", "top.fitTop1.M()", "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";m_{t}^{fit} [GeV], #eta^{j}<1.1; Permutations / 5 GeV" , 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("recoW1Mass_barrel" , "top.recoW1.M()" , "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";m_{W}^{reco} [GeV], #eta^{j}<1.1; Permutations / 5 GeV", 60,  0, 300));
    hists.push_back(MyHistogram("recoDeltaTopMass"    , "top.recoTop1.M()-top.recoTop2.M()"   , "", ";#Delta m_{t}^{reco} [GeV]; Permutations / 5 GeV", 80, -200, 200));
    hists.push_back(MyHistogram("recoDeltaTopMassBest", "top.recoTop1[0].M()-top.recoTop2[0].M()", "", ";#Delta m_{t}^{reco} [GeV]; Events / 5 GeV"      , 80, -200, 200));

    hists.push_back(MyHistogram("recoTopMass"        , std::vector<std::string>({"top.recoTop1.M()"   , "top.recoTop2.M()"   }), "", ";m_{t}^{reco} [GeV]; Top Quarks / 5 GeV", 70, 50, 400));
    hists.push_back(MyHistogram("recoTopMassBest"    , std::vector<std::string>({"top.recoTop1[0].M()", "top.recoTop2[0].M()"}), "", ";m_{t}^{reco} [GeV]; Top Quarks / 5 GeV", 70, 50, 400));
    hists.push_back(MyHistogram("recoTopMassPeak"    , std::vector<std::string>({"top.recoTop1.M()"   , "top.recoTop2.M()"   }), "", ";m_{t}^{reco} [GeV]; Top Quarks / 2 GeV", 45, 125, 215));
    hists.push_back(MyHistogram("recoTopMassPeakBest", std::vector<std::string>({"top.recoTop1[0].M()", "top.recoTop2[0].M()"}), "", ";m_{t}^{reco} [GeV]; Top Quarks / 2 GeV", 45, 125, 215));

    hists.push_back(MyHistogram("fitTop1Mass_vs_nVertex", "weight.nVertex", "top.fitTop1.M()", "", ";N_{Vertex}; m_{t}^{fit} [GeV]" , 8, 0, 40, 70, 50, 400));
    hists.push_back(MyHistogram("fitTop1Mass_vs_fitProb", "top.fitProb"   , "top.fitTop1.M()", "", ";P_{gof}; m_{t}^{fit} [GeV]"    , 10, 0,  1, 70, 50, 400));
    hists.push_back(MyHistogram("fitTop1Mass_vs_nJet", "jet.@jet.size()", "top.fitTop1.M()", "", ";N_{jet}; m_{t}^{fit} [GeV]", 11, 4, 15, 70, 50, 400));
    hists.push_back(MyHistogram("fitTop1Mass_vs_dRbb", "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))", "top.fitTop1.M()", "", ";#DeltaR_{b#bar{b}}; m_{t}^{fit} [GeV]" , 10, 1, 6, 70, 50, 400));
    hists.push_back(MyHistogram("recoW1Mass_vs_nVertex" , "weight.nVertex", "top.recoW1.M()" , "", ";N_{Vertex}; m_{W}^{reco} [GeV]", 8, 0, 40, 58, 10, 300));
    hists.push_back(MyHistogram("recoW1Mass_vs_fitProb" , "top.fitProb"   , "top.recoW1.M()" , "", ";P_{gof}; m_{W}^{reco} [GeV]"   , 10, 0,  1, 58, 10, 300));
    hists.push_back(MyHistogram("recoW1Mass_vs_nJet", "jet.@jet.size()", "top.recoW1.M()", "", ";N_{jet}; m_{W}^{reco} [GeV]", 11, 4, 15, 58, 10, 300));
    hists.push_back(MyHistogram("recoW1Mass_vs_jet3Pt", "jet.jet[3].Pt()", "top.recoW1.M()", "", ";p_{T}^{3} [GeV]; m_{W}^{reco} [GeV]", 7, 30, 100, 58, 10, 300));
    hists.push_back(MyHistogram("recoW1Mass_vs_dRqq", "sqrt(pow(top.recoW1Prod1.Eta()-top.recoW1Prod2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.recoW1Prod1.Phi()-top.recoW1Prod2.Phi()),2))", "top.recoW1.M()", "", ";#DeltaR_{q#bar{q}}; m_{W}^{reco}" , 10, 0, 5, 58, 10, 300));
    hists.push_back(MyHistogram("recoW1Mass_vs_dRbb", "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))", "top.recoW1.M()", "", ";#DeltaR_{b#bar{b}}; m_{W}^{reco}" , 10, 1, 6, 58, 10, 300));
  }
  if(plotSelectedForPlotting.find("WMassComponents")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("recoW1Prod1Mass", "top.recoW1Prod1.M()", "", ";m_{q} [GeV]; Permutations", 60, 0, 60));
    hists.push_back(MyHistogram("recoW1Prod2Mass", "top.recoW1Prod2.M()", "", ";m_{#bar{q}} [GeV]; Permutations", 60, 0, 60));
    hists.push_back(MyHistogram("recoW1SqrtProdEnergies", "sqrt(top.recoW1Prod1.E()*top.recoW1Prod2.E())", "", ";#sqrt{E_{q} E_{#bar{q}}} [GeV]; Permutations", 40, 0, 400));
    hists.push_back(MyHistogram("recoW1SqrtProdMomenta", "sqrt(top.recoW1Prod1.Px()*top.recoW1Prod2.Px() + top.recoW1Prod1.Py()*top.recoW1Prod2.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod2.Pz())", "", ";#sqrt{#vec{q} #vec{#bar{q}}} [GeV]; Permutations", 60, 0, 300));
    hists.push_back(MyHistogram("recoW1SqrtProdAbsMomenta", "sqrt(sqrt(top.recoW1Prod1.Px()*top.recoW1Prod1.Px() + top.recoW1Prod1.Py()*top.recoW1Prod1.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod1.Pz())*sqrt(top.recoW1Prod2.Px()*top.recoW1Prod2.Px() + top.recoW1Prod2.Py()*top.recoW1Prod2.Py() + top.recoW1Prod2.Pz()*top.recoW1Prod2.Pz()))", "", ";#sqrt{|#vec{q}| |#vec{#bar{q}}|} [GeV]; Permutations", 40, 0, 400));
    hists.push_back(MyHistogram("recoW1ProdAngle", "TMath::ACos((top.recoW1Prod1.Px()*top.recoW1Prod2.Px() + top.recoW1Prod1.Py()*top.recoW1Prod2.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod2.Pz())/sqrt((top.recoW1Prod1.Px()*top.recoW1Prod1.Px() + top.recoW1Prod1.Py()*top.recoW1Prod1.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod1.Pz())*(top.recoW1Prod2.Px()*top.recoW1Prod2.Px() + top.recoW1Prod2.Py()*top.recoW1Prod2.Py() + top.recoW1Prod2.Pz()*top.recoW1Prod2.Pz())))", "", ";#angle(q#bar{q}); Permutations", 32, 0, 3.2));
    hists.push_back(MyHistogram("recoW1MT", "top.recoW1.Mt()", "", ";m_{T,q} [GeV]; Permutations", 60, 0, 600));
    hists.push_back(MyHistogram("recoW1Prod1MT", "top.recoW1Prod1.Mt()", "", ";m_{T,q} [GeV]; Permutations", 60, 0, 300));
    hists.push_back(MyHistogram("recoW1Prod2MT", "top.recoW1Prod2.Mt()", "", ";m_{T,#bar{q}} [GeV]; Permutations", 60, 0, 300));
    hists.push_back(MyHistogram("recoW1ProdPtOverLeptonPt", "(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())/top.recoW2Prod1.Pt()", "", ";(p_{T,q}+p_{T,#bar{q}})/p_{T,lep}; Permutations", 50, 0, 10));
    hists.push_back(MyHistogram("recoW1Mass0", "sqrt(2*(sqrt(pow(top.recoW1Prod1.E(),2)-top.recoW1Prod1.M2())*sqrt(pow(top.recoW1Prod2.E(),2)-top.recoW1Prod2.M2()) - (top.recoW1Prod1.Px()*top.recoW1Prod2.Px() + top.recoW1Prod1.Py()*top.recoW1Prod2.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod2.Pz()) ))", "", ";m_{W,0}^{reco} [GeV]; Permutations", 60, 0, 300));
    hists.push_back(MyHistogram("recoW1Mass0Zoom", "sqrt(2*(sqrt(pow(top.recoW1Prod1.E(),2)-top.recoW1Prod1.M2())*sqrt(pow(top.recoW1Prod2.E(),2)-top.recoW1Prod2.M2()) - (top.recoW1Prod1.Px()*top.recoW1Prod2.Px() + top.recoW1Prod1.Py()*top.recoW1Prod2.Py() + top.recoW1Prod1.Pz()*top.recoW1Prod2.Pz()) ))", "", ";m_{W,0}^{reco} [GeV]; Permutations", 50, 60, 110));
  }
  // masses
  if(plotSelectedForPlotting.find("KinFitPulls")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("B1PullPt"   , "top.fitB1.Pt()-jet.jet[top.recoJetIdxB1].Pt()", "", ";#Delta p_{T}^{B1} [GeV]; Permutations / 5 GeV", 30, -30, 30));
    hists.push_back(MyHistogram("B2PullPt"   , "top.fitB2.Pt()-jet.jet[top.recoJetIdxB2].Pt()", "", ";#Delta p_{T}^{B2} [GeV]; Permutations / 5 GeV", 30, -30, 30));
    hists.push_back(MyHistogram("W1P1PullPt" , "top.fitW1Prod1.Pt()-jet.jet[top.recoJetIdxW1Prod1].Pt()", "", ";#Delta p_{T}^{W1L1} [GeV]; Permutations / 5 GeV", 20, -20, 20));
    hists.push_back(MyHistogram("W1P2PullPt" , "top.fitW1Prod2.Pt()-jet.jet[top.recoJetIdxW1Prod2].Pt()", "", ";#Delta p_{T}^{W1L2} [GeV]; Permutations / 5 GeV", 20, -20, 20));
    hists.push_back(MyHistogram("W2P1PullPt" , "top.fitW2Prod1.Pt()-jet.jet[top.recoJetIdxW2Prod1].Pt()", "", ";#Delta p_{T}^{W2L1} [GeV]; Permutations / 5 GeV", 20, -20, 20));
    hists.push_back(MyHistogram("W2P2PullPt" , "top.fitW2Prod2.Pt()-jet.jet[top.recoJetIdxW2Prod2].Pt()", "", ";#Delta p_{T}^{W2L2} [GeV]; Permutations / 5 GeV", 20, -20, 20));

    // needs a full constructor std::vector<std::string> as the cast from the initializer list only works for 4 or more entries ... WTF
    hists.push_back(MyHistogram("BJetPullPt" , std::vector<std::string>({"top.fitB1.Pt()-jet.jet[top.recoJetIdxB1].Pt()",
                                                                         "top.fitB2.Pt()-jet.jet[top.recoJetIdxB2].Pt()"}) , "", ";#Delta p_{T}^{B} [GeV]; B Jets / 5 GeV", 30, -30, 30));
    hists.push_back(MyHistogram("LJetPullPt" , {"top.fitW1Prod1.Pt()-jet.jet[top.recoJetIdxW1Prod1].Pt()",
                                                "top.fitW1Prod2.Pt()-jet.jet[top.recoJetIdxW1Prod2].Pt()",
                                                "top.fitW2Prod1.Pt()-jet.jet[top.recoJetIdxW2Prod1].Pt()",
                                                "top.fitW2Prod2.Pt()-jet.jet[top.recoJetIdxW2Prod2].Pt()"} , "", ";#Delta p_{T}^{L} [GeV]; Light Jets / 5 GeV", 20, -20, 20));

    hists.push_back(MyHistogram("B1PullEta"   , "top.fitB1.Eta()-jet.jet[top.recoJetIdxB1].Eta()", "", ";#Delta #eta^{B1}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("B2PullEta"   , "top.fitB2.Eta()-jet.jet[top.recoJetIdxB2].Eta()", "", ";#Delta #eta^{B2}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("W1P1PullEta" , "top.fitW1Prod1.Eta()-jet.jet[top.recoJetIdxW1Prod1].Eta()", "", ";#Delta #eta^{W1L1}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("W1P2PullEta" , "top.fitW1Prod2.Eta()-jet.jet[top.recoJetIdxW1Prod2].Eta()", "", ";#Delta #eta^{W1L2}; Permutations / 1 #times10^{-3}", 40, -0.02, 0.02));
    hists.push_back(MyHistogram("W2P1PullEta" , "top.fitW2Prod1.Eta()-jet.jet[top.recoJetIdxW2Prod1].Eta()", "", ";#Delta #eta^{W2L1}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("W2P2PullEta" , "top.fitW2Prod2.Eta()-jet.jet[top.recoJetIdxW2Prod2].Eta()", "", ";#Delta #eta^{W2L2}; Permutations / 1 #times10^{-3}", 40, -0.02, 0.02));

    // needs a full constructor std::vector<std::string> as the cast from the initializer list only works for 4 or more entries ... WTF
    hists.push_back(MyHistogram("BJetPullEta" , std::vector<std::string>({"top.fitB1.Eta()-jet.jet[top.recoJetIdxB1].Eta()",
                                                                          "top.fitB2.Eta()-jet.jet[top.recoJetIdxB2].Eta()"}) , "", ";#Delta #eta^{B} [GeV]; B Jets / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("LJetPullEta" , {"top.fitW1Prod1.Eta()-jet.jet[top.recoJetIdxW1Prod1].Eta()",
                                                 "top.fitW1Prod2.Eta()-jet.jet[top.recoJetIdxW1Prod2].Eta()",
                                                 "top.fitW2Prod1.Eta()-jet.jet[top.recoJetIdxW2Prod1].Eta()",
                                                 "top.fitW2Prod2.Eta()-jet.jet[top.recoJetIdxW2Prod2].Eta()"} , "", ";#Delta #eta^{L} [GeV]; Light Jets / 1 #times10^{-3}", 40, -0.02, 0.02));

    hists.push_back(MyHistogram("B1PullPhi"   , "top.fitB1.Phi()-jet.jet[top.recoJetIdxB1].Phi()", "", ";#Delta #phi^{B1}; Permutations / 2.5 #times10^{-4}", 40, -0.005, 0.005));
    hists.push_back(MyHistogram("B2PullPhi"   , "top.fitB2.Phi()-jet.jet[top.recoJetIdxB2].Phi()", "", ";#Delta #phi^{B2}; Permutations / 2.5 #times10^{-4}", 40, -0.005, 0.005));
    hists.push_back(MyHistogram("W1P1PullPhi" , "top.fitW1Prod1.Phi()-jet.jet[top.recoJetIdxW1Prod1].Phi()", "", ";#Delta #phi^{W1L1}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("W1P2PullPhi" , "top.fitW1Prod2.Phi()-jet.jet[top.recoJetIdxW1Prod2].Phi()", "", ";#Delta #phi^{W1L2}; Permutations / 1 #times10^{-3}", 40, -0.02, 0.02));
    hists.push_back(MyHistogram("W2P1PullPhi" , "top.fitW2Prod1.Phi()-jet.jet[top.recoJetIdxW2Prod1].Phi()", "", ";#Delta #phi^{W2L1}; Permutations / 5 #times10^{-4}", 40, -0.01, 0.01));
    hists.push_back(MyHistogram("W2P2PullPhi" , "top.fitW2Prod2.Phi()-jet.jet[top.recoJetIdxW2Prod2].Phi()", "", ";#Delta #phi^{W2L2}; Permutations / 1 #times10^{-3}", 40, -0.02, 0.02));

    // needs a full constructor std::vector<std::string> as the cast from the initializer list only works for 4 or more entries ... WTF
    hists.push_back(MyHistogram("BJetPullPhi" , std::vector<std::string>({"top.fitB1.Phi()-jet.jet[top.recoJetIdxB1].Phi()",
                                                                          "top.fitB2.Phi()-jet.jet[top.recoJetIdxB2].Phi()"}) , "", ";#Delta #phi^{B} [GeV]; B Jets / 2.5 #times10^{-4}", 40, -0.005, 0.005));
    hists.push_back(MyHistogram("LJetPullPhi" , {"top.fitW1Prod1.Phi()-jet.jet[top.recoJetIdxW1Prod1].Phi()",
                                                 "top.fitW1Prod2.Phi()-jet.jet[top.recoJetIdxW1Prod2].Phi()",
                                                 "top.fitW2Prod1.Phi()-jet.jet[top.recoJetIdxW2Prod1].Phi()",
                                                 "top.fitW2Prod2.Phi()-jet.jet[top.recoJetIdxW2Prod2].Phi()"} , "", ";#Delta #phi^{L} [GeV]; Light Jets / 1 #times10^{-3}", 40, -0.02, 0.02));
  }
  // BReg cross check masses
  if(plotSelectedForPlotting.find("BRegCrossCheckExtraMasses")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTop1MassHighRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()>70", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassMediumRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()>50&&top.fitB1.Pt()<70", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassLowRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()<50", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassHighPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]>70", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassMediumPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]>50&&BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]<70", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassLowPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]<50", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();


    hists.push_back(MyHistogram("fitTop1MassCPHighRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()>70&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassCPMediumRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()>50&&top.fitB1.Pt()<70&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassCPLowRecoB"    , "top.fitTop1.M()"   , "top.fitB1.Pt()<50&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassCPHighPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]>70&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassCPMediumPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]>50&&BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]<70&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("fitTop1MassCPLowPtGenB"    , "top.fitTop1.M()"   , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]<50&&"+topBranchName_+"combinationType==1", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();


    hists.push_back(MyHistogram("B1GenJetPtCP"    , "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]",topBranchName_+"combinationType==1", ";p_{T}^{gen} [GeV]; Permutations / 5 GeV", 20, 0, 200));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("B1GenPartonPtCP"    , "BRegJet.genPartonPt["+topBranchName_+"recoJetIdxB1]",topBranchName_+"combinationType==1", ";p_{T}^{gen} [GeV]; Permutations / 5 GeV", 20, 0, 200));
    hists.back().SetFitGaussToCore();

    hists.push_back(MyHistogram("B1JetFitPtCP"    , "top.fitB1.Pt()",topBranchName_+"combinationType==1", ";p_{T}^{fit} [GeV]; Permutations / 5 GeV", 20, 0, 200));
    hists.back().SetFitGaussToCore();


    hists.push_back(MyHistogram("B1TopMassCPVsGenJetPt","BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]", topBranchName_+"fitTop1.M()", topBranchName_+"combinationType==1",";p_{T}^{gen}; m_{t}^{fit} [GeV]",20,0,200, 20, 50, 400 ));
    hists.back().ConfigurePlotRanges(0.95,1.05, 150, 200);
    hists.push_back(MyHistogram("B1TopMassCPVsGenPartonPt","BRegJet.genPartonPt["+topBranchName_+"recoJetIdxB1]", topBranchName_+"fitTop1.M()", topBranchName_+"combinationType==1",";p_{T}^{gen}; m_{t}^{fit} [GeV]",20,0,200, 20, 50, 400 ));
    hists.back().ConfigurePlotRanges(0.95,1.05, 150, 200);
    hists.push_back(MyHistogram("B1TopMassCPVsFitPt","top.fitB1.Pt()", topBranchName_+"fitTop1.M()", topBranchName_+"combinationType==1",";p_{T}^{fit}; m_{t}^{fit} [GeV]",20,0,200, 20, 50, 400 ));
    hists.back().ConfigurePlotRanges(0.95,1.05, 150, 200);
  }
  // light pulls
  if(plotSelectedForPlotting.find("LightPulls")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("pullW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pull[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{q_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod2_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),-(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta())))))", "", ";#theta_{pull,ch}^{q_{2} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{q_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));

    hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_0",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),2)) < 1.",
				";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_1",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),2)) > 1.",
				";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_W1Prod2_deltaR_2",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),2)) > 2.",
				";#theta_{pull,ch}^{q_{1} #rightarrow q_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));
  }

  // b pulls
  if(plotSelectedForPlotting.find("BJetPulls")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("pullChargedB1_B2", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))", "", ";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{b_{1} #rightarrow beam}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB2_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitB2.Phi()),-(top.fitB1.Eta()-top.fitB2.Eta())))))", "", ";#theta_{pull,ch}^{b_{2} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB2_beam", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB2].Phi()-(TMath::Pi()+TMath::ATan2(-1,-0))))", "", ";#theta_{pull,ch}^{b_{2} #rightarrow beam}; Permutations", 20, 0, 3.1416));

    hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_0",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(TVector2::Phi_mpi_pi(top.fitB2.Phi()-top.fitB1.Phi()),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) < 1.",
				";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_1",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(TVector2::Phi_mpi_pi(top.fitB2.Phi()-top.fitB1.Phi()),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) > 1.",
				";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_B2_deltaR_2",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(TVector2::Phi_mpi_pi(top.fitB2.Phi()-top.fitB1.Phi()),2)+pow(top.fitB2.Eta()-top.fitB1.Eta(),2)) > 2.",
				";#theta_{pull,ch}^{b_{1} #rightarrow b_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));
  }

  // mixed pulls
  if(plotSelectedForPlotting.find("MixedPulls")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("pullChargedW1Prod1_B1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))", "", ";#theta_{pull}^{q_{1} #rightarrow b_{1}}; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_W1Prod1", "abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitB1.Phi()),-(top.fitW1Prod1.Eta()-top.fitB1.Eta())))))", "", ";#theta_{pull}^{b_{1} #rightarrow q_{1}}; Permutations", 20, 0, 3.1416));

    hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_0",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod1.Phi()),2)) < 1.",
				";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_1",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod1.Phi()),2)) > 1.",
				";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedW1Prod1_B1_deltaR_2",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB1.Phi()-top.fitW1Prod1.Phi()),-(top.fitB1.Eta()-top.fitW1Prod1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod1.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod1.Phi()),2)) > 2.",
				";#theta_{pull,ch}^{q_{1} #rightarrow b_{1}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));

    hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_0",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod2.Phi()),2)) < 1.",
				";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR < 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_1",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod2.Phi()),2)) > 2.",
				";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR > 1; Permutations", 20, 0, 3.1416));
    hists.push_back(MyHistogram("pullChargedB1_W1Prod2_deltaR_2",
				"abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitB1.Phi()),-(top.fitW1Prod2.Eta()-top.fitB1.Eta())))))",
				"sqrt(pow(top.fitB1.Eta()-top.fitW1Prod2.Eta(),2)+pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitW1Prod2.Phi()),2)) > 2.",
				";#theta_{pull,ch}^{b_{1} #rightarrow q_{2}}, #DeltaR > 2; Permutations", 20, 0, 3.1416));
  }

  // Basic kinematics
  if(plotSelectedForPlotting.find("JetPts")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("jet1Pt", "jet.jet[0].Pt()", "", ";p_{T}^{jet 1} [GeV]; Events / 10 GeV", 45, 0, 450));
    hists.push_back(MyHistogram("jet2Pt", "jet.jet[1].Pt()", "", ";p_{T}^{jet 2} [GeV]; Events / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("jet3Pt", "jet.jet[2].Pt()", "", ";p_{T}^{jet 3} [GeV]; Events / 5 GeV", 50, 0, 250));
    hists.push_back(MyHistogram("jet4Pt", "jet.jet[3].Pt()", "", ";p_{T}^{jet 4} [GeV]; Events / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("jet5Pt", "jet.jet[4].Pt()", "jet.jet[4].Pt()>0", ";p_{T}^{jet 5} [GeV]; Events / 3 GeV", 50, 0, 150));
    hists.push_back(MyHistogram("jet6Pt", "jet.jet[5].Pt()", "jet.jet[5].Pt()>0", ";p_{T}^{jet 6} [GeV]; Events / 2 GeV", 50, 0, 100));
    hists.push_back(MyHistogram("lJetPt", "jet.jet.Pt()", "jet.bTagCSV<0.679", ";p_{T} [GeV], untagged (CSVM); Jets / 10 GeV", 45, 0, 450));
    hists.push_back(MyHistogram("bJetPt", "jet.jet.Pt()", "jet.bTagCSV>0.679", ";p_{T} [GeV], b-tagged (CSVM); Jets / 10 GeV", 45, 0, 450));
    hists.push_back(MyHistogram("lJetPt_vs_nVertex" , "weight.nVertex", "jet.jet.Pt()" , "jet.bTagCSV<0.679", ";N_{Vertex}; p_{T} [GeV], untagged (CSVM)", 8, 0, 40, 44, 10, 450));
    hists.push_back(MyHistogram("bJetPt_vs_nVertex" , "weight.nVertex", "jet.jet.Pt()" , "jet.bTagCSV>0.679", ";N_{Vertex}; p_{T} [GeV], b-tagged (CSVM)", 8, 0, 40, 44, 10, 450));
  }
  
  if(plotSelectedForPlotting.find("JetEtas")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("jet1Eta", "jet.jet[0].Eta()", "", ";#eta^{1}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet2Eta", "jet.jet[1].Eta()", "", ";#eta^{2}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet3Eta", "jet.jet[2].Eta()", "", ";#eta^{3}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet4Eta", "jet.jet[3].Eta()", "", ";#eta^{4}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet5Eta", "jet.jet[4].Eta()", "jet.jet[4].Pt()>0", ";#eta^{5}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet6Eta", "jet.jet[5].Eta()", "jet.jet[5].Pt()>0", ";#eta^{6}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("lJetEta", "jet.jet.Eta()", "jet.bTagCSV<0.679", ";#eta, untagged (CSVM); Jets / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("bJetEta", "jet.jet.Eta()", "jet.bTagCSV>0.679", ";#eta, b-tagged (CSVM); Jets / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("lJetEta_vs_nVertex" , "weight.nVertex", "abs(jet.jet.Eta())" , "jet.bTagCSV<0.679", ";N_{Vertex}; #eta, untagged (CSVM)", 6, 0, 30, 30, 0, 3));
    hists.push_back(MyHistogram("bJetEta_vs_nVertex" , "weight.nVertex", "abs(jet.jet.Eta())" , "jet.bTagCSV>0.679", ";N_{Vertex}; #eta, b-tagged (CSVM)", 6, 0, 30, 30, 0, 3));
  }

  if(plotSelectedForPlotting.find("JetPhis")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("jet1Phi", "jet.jet[0].Phi()", "", ";#phi^{1}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet2Phi", "jet.jet[1].Phi()", "", ";#phi^{2}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet3Phi", "jet.jet[2].Phi()", "", ";#phi^{3}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet4Phi", "jet.jet[3].Phi()", "", ";#phi^{4}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet5Phi", "jet.jet[4].Phi()", "jet.jet[4].Pt()>0", ";#phi^{5}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("jet6Phi", "jet.jet[5].Phi()", "jet.jet[5].Pt()>0", ";#phi^{6}; Events / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("lJetPhi", "jet.jet.Phi()", "jet.bTagCSV<0.679", ";#phi, untagged (CSVM); Jets / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("bJetPhi", "jet.jet.Phi()", "jet.bTagCSV>0.679", ";#phi, b-tagged (CSVM); Jets / 0.2", 30, -3, 3));
    hists.push_back(MyHistogram("lJetPhi_vs_nVertex" , "weight.nVertex", "abs(jet.jet.Phi())" , "jet.bTagCSV<0.679", ";N_{Vertex}; #phi, untagged (CSVM)", 6, 0, 30, 30, 0, 3));
    hists.push_back(MyHistogram("bJetPhi_vs_nVertex" , "weight.nVertex", "abs(jet.jet.Phi())" , "jet.bTagCSV>0.679", ";N_{Vertex}; #phi, b-tagged (CSVM)", 6, 0, 30, 30, 0, 3));
  }
  
  if(plotSelectedForPlotting.find("JetMass")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("jet1Mass", "jet.jet[0].M()", "", ";m^{1}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("jet2Mass", "jet.jet[1].M()", "", ";m^{2}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("jet3Mass", "jet.jet[2].M()", "", ";m^{3}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("jet4Mass", "jet.jet[3].M()", "", ";m^{4}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("jet5Mass", "jet.jet[4].M()", "jet.jet[4].Pt()>0", ";m^{5}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("jet6Mass", "jet.jet[5].M()", "jet.jet[5].Pt()>0", ";m^{6}; Events / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("lJetMass", "jet.jet.M()", "jet.bTagCSV<0.679", ";m, untagged (CSVM); Jets / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("bJetMass", "jet.jet.M()", "jet.bTagCSV>0.679", ";m, b-tagged (CSVM); Jets / 0.2", 50, 0, 50));
    hists.push_back(MyHistogram("lJetMass_vs_nVertex" , "weight.nVertex", "abs(jet.jet.M())" , "jet.bTagCSV<0.679", ";N_{Vertex}; m, untagged (CSVM)", 6, 0, 30, 50, 0, 50));
    hists.push_back(MyHistogram("bJetMass_vs_nVertex" , "weight.nVertex", "abs(jet.jet.M())" , "jet.bTagCSV>0.679", ";N_{Vertex}; m, b-tagged (CSVM)", 6, 0, 30, 50, 0, 50));
  }

  if(plotSelectedForPlotting.find("TopPts")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTop1Pt"    , "top.fitTop1.Pt()"   , "", ";p_{T,t}^{fit} [GeV]; Permutations / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("fitTop1PtBest", "top.fitTop1[0].Pt()", "", ";p_{T,t}^{fit} [GeV]; Events / 10 GeV"      , 40, 0, 400));
    hists.push_back(MyHistogram("fitTop1Pt_barrel", "top.fitTop1.Pt()", "top.fitB1.Eta()<1.1 & top.fitW1Prod1.Eta()<1.1 & top.fitW1Prod2.Eta()<1.1", ";p_{T,t}^{fit} [GeV], #eta^{j}<1.1; Permutations / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("fitRlb", "(top.fitB1.Pt()+top.fitB2.Pt())/(top.fitW1Prod1.Pt()+top.fitW1Prod2.Pt())", "", ";R_{lb}^{fit}; Permutations", 40, 0, 4));
    hists.push_back(MyHistogram("recoRlb", "(top.recoB1.Pt()+top.recoB2.Pt())/(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())", "", ";R_{lb}^{reco}; Permutations", 40, 0, 4));
    hists.push_back(MyHistogram("fitTop1Pt_vs_fitTop1Mass"    , "top.fitTop1.Pt()"   , "top.fitTop1.M()"   , "", ";p_{T,t}^{fit} [GeV];m_{t}^{fit} [GeV]", 16, 0, 400, 70, 50, 400));
    hists.push_back(MyHistogram("fitTop1Pt_vs_fitTop1MassBest", "top.fitTop1[0].Pt()", "top.fitTop1[0].M()", "", ";p_{T,t}^{fit} [GeV];m_{t}^{fit} [GeV]", 16, 0, 400, 70, 50, 400));
    hists.push_back(MyHistogram("fitTopPt"    , std::vector<std::string>({"top.fitTop1.Pt()"   , "top.fitTop2.Pt()"   }) , "", ";p_{T,t}^{fit} [GeV]; Top Quarks / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("fitTopPtBest", std::vector<std::string>({"top.fitTop1[0].Pt()", "top.fitTop2[0].Pt()"}) , "", ";p_{T,t}^{fit} [GeV]; Top Quarks / 10 GeV", 40, 0, 400));
  }

  if(plotSelectedForPlotting.find("TTBarSystem")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTTBarPt"    , "top.fitTTBar.Pt()"          , "", ";p_{T,t#bar{t}}^{fit} [GeV]; Permutations / 10 GeV", 40, 0,  400));
    hists.push_back(MyHistogram("fitTTBarPtBest", "top.fitTTBar[0].Pt()"       , "", ";p_{T,t#bar{t}}^{fit} [GeV]; Events / 10 GeV"      , 40, 0,  400));
    hists.push_back(MyHistogram("fitTTBarM"     , "top.fitTTBar.M()"           , "", ";m_{t#bar{t}}^{fit} [GeV]; Permutations / 30 GeV"  , 40, 300, 1500));
    hists.push_back(MyHistogram("fitTTBarMBest" , "top.fitTTBar[0].M()"        , "", ";m_{t#bar{t}}^{fit} [GeV]; Events / 30 GeV"        , 40, 300, 1500));
    hists.push_back(MyHistogram("fitTTBarY"     , "top.fitTTBar.Rapidity()"    , "", ";y_{t#bar{t}}^{fit}; Permutations / 0.2"           , 40, -4, 4));
    hists.push_back(MyHistogram("fitTTBarYBest" , "top.fitTTBar[0].Rapidity()" , "", ";y_{t#bar{t}}^{fit}; Events / 0.2"                 , 40, -4, 4));
  }

  if(plotSelectedForPlotting.find("newCuts")!=plotSelectedForPlotting.end()){
	    hists.push_back(MyHistogram("recoHt"    , "top.recoB1.Pt()+top.recoB2.Pt()+top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt()+top.recoW2Prod1.Pt()+top.recoW2Prod2.Pt()+jet.met.Pt()" , "", ";H_{T}^{reco} [GeV]; Permutations / 100 GeV", 30, 0,  3000));
	    hists.push_back(MyHistogram("recoHtBest", "top.recoB1[0].Pt()+top.recoB2[0].Pt()+top.recoW1Prod1[0].Pt()+top.recoW1Prod2[0].Pt()+top.recoW2Prod1[0].Pt()+top.recoW2Prod2[0].Pt()+jet.met.Pt()" , "", ";H_{T}^{reco} [GeV]; Events / 100 GeV", 30, 0,  3000));
	    hists.push_back(MyHistogram("recoDeltaRLepToB"    , "sqrt(pow(top.recoB1.Eta()-top.recoW2Prod1.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.recoB1.Phi()-top.recoW2Prod1.Phi()),2))" , "", ";/{DelraR}_{/{mu} ,b_{lep}}^{reco} ; Permutations / 0.1", 30, 0,  3));
	    hists.push_back(MyHistogram("recoDeltaRLepToBBest"    , "sqrt(pow(top.recoB1[0].Eta()-top.recoW2Prod1[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.recoB1[0].Phi()-top.recoW2Prod1[0].Phi()),2))" , "", ";/{DelraR}_{/{mu} ,b_{lep}}^{reco} ; Events / 0.1", 30, 0,  3));

  }

  if(plotSelectedForPlotting.find("WBosonExtra")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitWPt"     , std::vector<std::string>({"top.fitW1.Pt()"         , "top.fitW2.Pt()"         }) , "", ";p_{T,W}^{fit} [GeV]; W Bosons / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("fitWPtBest" , std::vector<std::string>({"top.fitW1[0].Pt()"      , "top.fitW2[0].Pt()"      }) , "", ";p_{T,W}^{fit} [GeV]; W Bosons / 10 GeV", 40, 0, 400));
    hists.push_back(MyHistogram("fitWY"      , std::vector<std::string>({"top.fitW1.Rapidity()"   , "top.fitW2.Rapidity()"   }) , "", ";y_{W}^{fit}; W Bosons / 0.2"           , 40, -4, 4));
    hists.push_back(MyHistogram("fitWYBest"  , std::vector<std::string>({"top.fitW1[0].Rapidity()", "top.fitW2[0].Rapidity()"}) , "", ";y_{W}^{fit}; W Bosons / 0.2"           , 40, -4, 4));
  }

  if(plotSelectedForPlotting.find("TopExtra")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitTopY"      , std::vector<std::string>({"top.fitTop1.Rapidity()"   , "top.fitTop2.Rapidity()"   }) , "", ";y_{t}^{fit}; Top Quarks / 0.2", 40, -4, 4));
    hists.push_back(MyHistogram("fitTopYBest"  , std::vector<std::string>({"top.fitTop1[0].Rapidity()", "top.fitTop2[0].Rapidity()"}) , "", ";y_{t}^{fit}; Top Quarks / 0.2", 40, -4, 4));
  }


  if(plotSelectedForPlotting.find("FitPts")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitB1Pt", "top.fitB1.Pt()", "", ";p_{T}^{B1} [GeV]; Permutations / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("fitB2Pt", "top.fitB2.Pt()", "", ";p_{T}^{B2} [GeV]; Permutations / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("fitW1Prod1Pt", "top.fitW1Prod1.Pt()", "", ";p_{T}^{W1,1} [GeV]; Permutations / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("fitW1Prod2Pt", "top.fitW1Prod2.Pt()", "", ";p_{T}^{W1,2} [GeV]; Permutations / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("fitW2Prod1Pt", "top.fitW2Prod1.Pt()", "", ";p_{T}^{W2,1} [GeV]; Permutations / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("fitW2Prod2Pt", "top.fitW2Prod2.Pt()", "", ";p_{T}^{W2,2} [GeV]; Permutations / 5 GeV", 40, 0, 200));
  }

  // lepton+jets
  if(plotSelectedForPlotting.find("LeptonJets")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("leptonPt", "top.recoW2Prod1[0].Pt()", "", ";p_{T}^{lepton} [GeV]; Events / 5 GeV", 40, 0, 200));
    hists.push_back(MyHistogram("leptonEta", "top.recoW2Prod1[0].Eta()", "", ";#eta^{lepton}; Events", 25, -2.5, 2.5));
    hists.push_back(MyHistogram("leptonPhi", "top.recoW2Prod1[0].Phi()", "", ";#phi^{lepton}; Events", 30, -3, 3));
    hists.push_back(MyHistogram("leptonDRJet", (TString("Min$(")+MyMa::deltaR("top.recoW2Prod1[0]", "jet.jet")+TString(")")).Data(), "", ";#DeltaR(l,nearest jet); Events", 40, 0, 4));
    hists.push_back(MyHistogram("MET", "top.recoW2Prod2[0].Pt()", "", ";E_{T}^{miss} [GeV]; Events / 5 GeV", 40, 0, 200));
  }
  
  if(plotSelectedForPlotting.find("LeptonJetsExtra")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("fitW1Pt", "top.fitW1.Pt()"   , "", ";p_{T,W,had}^{fit} [GeV]; Permutations", 40, 0, 400));
    hists.push_back(MyHistogram("fitW2Pt", "top.fitW2.Pt()"   , "", ";p_{T,W,lep}^{fit} [GeV]; Permutations", 40, 0, 400));
    hists.push_back(MyHistogram("fitTop1MassAtlas"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations / 5 GeV", 90, 130, 220));
    hists.push_back(MyHistogram("recoW1MassZoom"    , "top.recoW1.M()"   , "", ";m_{W}^{reco} [GeV]; Permutations / 2 GeV", 50, 60, 110));
    //	  hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("recoW2Mass"    , "top.recoW2.M()"   , "", ";m_{W,lep}^{reco} [GeV]; Permutations / 5 GeV", 60, 0, 300));
  }
  
    if(plotSelectedForPlotting.find("LeptonJetsSearch")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("SearchFitTop1Mass"    , "top.fitTop1.M()"   , "", ";m_{t}^{fit} [GeV]; Permutations / 20 GeV", 50, 300, 1300));
    hists.push_back(MyHistogram("SearchMET", "top.recoW2Prod2[0].Pt()", "", ";p_{T}^{miss} [GeV]; Events", 50, 0, 1000));
    hists.push_back(MyHistogram("SearchNu", "top.fitW2Prod2[0].Pt()", "", ";p_{T}^{#nu} [GeV]; Events", 50, 0, 1000));
  }
  
  if(plotSelectedForPlotting.find("ChargedHiggs")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("recoW1MassTopCut"    , "top.recoW1.M()"   , "abs(top.recoTop1.M()-172.5) < 10.", ";m_{W}^{reco} [GeV]; Permutations / 2 GeV", 100, 0, 200));
  }

  // jet details
  if(plotSelectedForPlotting.find("JetDetails")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("nChargedHadronsBJetT", "jet.nChargedHadrons", "jet.bTagCSV>0.898", ";N_{ch}; B Jets", 40, 0, 200));
    hists.push_back(MyHistogram("nChargedHadronsBJetM", "jet.nChargedHadrons", "jet.bTagCSV>0.679", ";N_{ch}; B Jets", 40, 0, 200));
    hists.push_back(MyHistogram("nChargedHadronsLJetT", "jet.nChargedHadrons", "jet.bTagCSV<0.898", ";N_{ch}; Non B Jets", 40, 0, 200));
    hists.push_back(MyHistogram("nChargedHadronsLJetM", "jet.nChargedHadrons", "jet.bTagCSV<0.679", ";N_{ch}; Non B Jets", 40, 0, 200));
  }
  
  // b-tagging
  if(plotSelectedForPlotting.find("JetBTag")!=plotSelectedForPlotting.end()){
	  hists.push_back(MyHistogram("nBJet4", "(jet.bTagCSV[0]>0.679) + (jet.bTagCSV[1]>0.679) + (jet.bTagCSV[2]>0.679) + (jet.bTagCSV[3]>0.679)", "", ";N_{b jet}^{4} (CSVM); Events", 5, -0.5, 4.5));
	  hists.push_back(MyHistogram("nBJet", "Sum$((jet.jet.Pt()>30 & jet.bTagCSV>0.679))", "", ";N_{b jet} (CSVM); Events", 5, -0.5, 4.5));
	  hists.push_back(MyHistogram("bTagCSV", "jet.bTagCSV", "jet.jet.Pt()>30", ";CSV discriminator; Jets", 41, -1, 1.05));
  }

  // event observables
  if(plotSelectedForPlotting.find("EventObservables")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("nJet", "jet.@jet.size()", "", ";N_{jet}; Events", 15, -0.5, 14.5));
    hists.push_back(MyHistogram("nJet30", "Sum$((jet.jet.Pt()>30))", "", ";N_{jet}; Events", 15, -0.5, 14.5));
    hists.push_back(MyHistogram("nVertex", "weight.nVertex", "", ";N_{PV}; Events", 35, 0, 35));
  }

  if(plotSelectedForPlotting.find("eventTreeWithJEC")!=plotSelectedForPlotting.end()){//50,60,110
    hists.push_back(MyHistogram("RecoW1MassTreeCorr",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "1.", "1.").Data(),"","; m_{W}^{treeCorr}; Permutations / 2.5 GeV",28,50,120));
    std::cout << MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "1.", "1.").Data() << std::endl;
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassStandardL1",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{StandardL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassNoPtL1",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{noPtL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassNoPtL1PtCut",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod1]*top.recoW1Prod1.Pt()>30&&JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod2]*top.recoW1Prod2.Pt()>30","; m_{W}^{noPtL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassPtFix",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{pt dep. L3Residual}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassPtFixV2",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "JetCorrL1L2L3PlResPtFixV2[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResPtFixV2[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{pt dep. L3Residual V2}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();

    hists.push_back(MyHistogram("RecoTop1MassPtFix",MyMa::invariantMassThreeVec("top.recoW1Prod1", "top.recoW1Prod2", "top.recoB1", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxB1]/TreeCorrFactor[top.recoJetIdxB1]").Data(),"","; m_{top, reco}^{pt dep. L3Residual}; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoTop1MassPtFixV2",MyMa::invariantMassThreeVec("top.recoW1Prod1", "top.recoW1Prod2", "top.recoB1", "JetCorrL1L2L3PlResPtFixV2[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlResPtFixV2[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]", "JetCorrL1L2L3PlResPtFixV2[top.recoJetIdxB1]/TreeCorrFactor[top.recoJetIdxB1]").Data(),"","; m_{top, reco}^{pt dep. L3Residual V2}; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoTop1MassStandardL1",MyMa::invariantMassThreeVec("top.recoW1Prod1", "top.recoW1Prod2", "top.recoB1", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod1]/TreeCorrFactor[top.recoJetIdxW1Prod1]", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod2]/TreeCorrFactor[top.recoJetIdxW1Prod2]", "JetCorrL1L2L3PlRes[top.recoJetIdxB1]/TreeCorrFactor[top.recoJetIdxB1]").Data(),"","; m_{top, reco}^{StandardL1}; Permutations / 5 GeV", 70, 50, 400));
    hists.back().SetFitGaussToCore();



    hists.push_back(MyHistogram("RecoW1MassBestStandardL1",MyMa::invariantMass("top.recoW1Prod1[0]", "top.recoW1Prod2[0]", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod1[0]]/TreeCorrFactor[top.recoJetIdxW1Prod1[0]]", "JetCorrL1L2L3PlRes[top.recoJetIdxW1Prod2[0]]/TreeCorrFactor[top.recoJetIdxW1Prod2[0]]").Data(),"","; m_{W}^{StandardL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassBestNoPtL1",MyMa::invariantMass("top.recoW1Prod1[0]", "top.recoW1Prod2[0]", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod1[0]]/TreeCorrFactor[top.recoJetIdxW1Prod1[0]]", "JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod2[0]]/TreeCorrFactor[top.recoJetIdxW1Prod2[0]]").Data(),"","; m_{W}^{noPtL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassBestPtFix",MyMa::invariantMass("top.recoW1Prod1[0]", "top.recoW1Prod2[0]", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod1[0]]/TreeCorrFactor[top.recoJetIdxW1Prod1[0]]", "JetCorrL1L2L3PlResPtFix[top.recoJetIdxW1Prod2[0]]/TreeCorrFactor[top.recoJetIdxW1Prod2[0]]").Data(),"","; m_{W}^{pt dep. L3Residual}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();



    hists.push_back(MyHistogram("RecoW1MassUseMCL1",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "(JetCorrL2L3L1MC[top.recoJetIdxW1Prod1]*JetCorrL2L3Res[top.recoJetIdxW1Prod1])/TreeCorrFactor[top.recoJetIdxW1Prod1]", "(JetCorrL2L3L1MC[top.recoJetIdxW1Prod2]*JetCorrL2L3Res[top.recoJetIdxW1Prod2])/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{MCL1}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1MassUseMCL1NoResidual",MyMa::invariantMass("top.recoW1Prod1", "top.recoW1Prod2", "(JetCorrL2L3L1MC[top.recoJetIdxW1Prod1])/TreeCorrFactor[top.recoJetIdxW1Prod1]", "(JetCorrL2L3L1MC[top.recoJetIdxW1Prod2])/TreeCorrFactor[top.recoJetIdxW1Prod2]").Data(),"","; m_{W}^{MCL1noResidual}; Permutations / 2.5 GeV",28,50,120));
    hists.back().SetFitGaussToCore();


    hists.push_back(MyHistogram("RecoW1Prod1Pt","top.recoW1Prod1.Pt()","","; p_{T,W}^{treeCorr}; Permutations / 2.5 GeV",50,20,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1Prod1PtNoPtL1","JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod1]*top.recoW1Prod1.Pt()","","; p_{T,W}^{noPtL1}; Permutations / 2.5 GeV",50,20,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1Prod2Pt","top.recoW1Prod2.Pt()","","; p_{T,W}^{treeCorr}; Permutations / 2.5 GeV",50,20,120));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoW1Prod2PtNoPtL1","JetCorrL1L2L3PlResNoPt[top.recoJetIdxW1Prod2]*top.recoW1Prod2.Pt()","","; p_{T,W}^{noPtL1}; Permutations / 2.5 GeV",50,20,120));
    hists.back().SetFitGaussToCore();

  }
  if(plotSelectedForPlotting.find("BRegExtraPlots")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo; 
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("BothBs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),
				  std::vector<std::string>({HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"),
					HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB2]")}),
				  "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));

      hists.push_back(MyHistogram("BothBsLogY"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),
				  std::vector<std::string>({HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"),
					HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB2]")}),
				  "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().ConfigureExtraOptions(false, "1.", false, "", false, true);


      hists.push_back(MyHistogram("B1LogY"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().ConfigureExtraOptions(false, "1.", false, "", false, true);


      hists.push_back(MyHistogram("B1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      //		  hists.back().SetFitGaussToCore();

      //                                                "top.fitW2Prod2.Phi()-jet.jet[top.recoJetIdxW2Prod2].Phi()"}

      //		  hists.back().SetFitGaussToCore();
      //		  hists.push_back(MyHistogram("CorrFactorWeightedB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      //		  hists.back().ConfigureExtraOptions(false, "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]");
      //		  hists.push_back(MyHistogram("GBRCorrFactorWeightedB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      //		  hists.back().ConfigureExtraOptions(false, "BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1]");
      //		  hists.push_back(MyHistogram("B1HighGBR"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1]>1.6",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      //		  hists.push_back(MyHistogram("W1Prod1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxW1Prod1]"), "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      //		  //2D
      hists.push_back(MyHistogram("B1TopMassVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), topBranchName_+"fitTop1.M()", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; m_{t}^{fit} [GeV]",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 50, 400 ));
      hists.back().ConfigurePlotRanges(0.95,1.05, 150, 250);
      hists.push_back(MyHistogram("B1TopMassCPVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), topBranchName_+"fitTop1.M()", topBranchName_+"combinationType==1",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; m_{t}^{fit} [GeV]",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 50, 400 ));
      hists.back().ConfigurePlotRanges(0.95,1.05, 150, 200);

      //		  hists.push_back(MyHistogram("B1GBRTrainFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1]", "",";"+helperMyBRegVarInfo.varNames.at(bvar_i)+"; GBRFactor",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0., 2. ));
      //		  hists.push_back(MyHistogram("B1recoRlbVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "(top.recoB1.Pt()+top.recoB2.Pt())/(top.recoW1Prod1.Pt()+top.recoW1Prod2.Pt())", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";R_{lb}^{reco}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 4));
      hists.push_back(MyHistogram("B1TrueCorrFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{gen}/p_{T}^{corr}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
      hists.push_back(MyHistogram("B1OnlyCPTrueCorrFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]", "top.combinationType==1", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{gen}/p_{T}^{corr}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
      hists.push_back(MyHistogram("B1TrueCorrFactorGBRCorrectedVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/(BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]*BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1])", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{gen}/p_{T}^{corr}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
      hists.push_back(MyHistogram("B1OnlyCPTrueCorrFactorGBRCorrectedVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/(BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]*BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1])", "top.combinationType==1", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{gen}/p_{T}^{corr}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
      hists.push_back(MyHistogram("B1FitTrueCorrFactorVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),HelperFunctions::addProperArrayIndex(helperMyBRegVarInfo.varForms.at(bvar_i),"["+topBranchName_+"recoJetIdxB1]"), "bRegTop.fitB1.Pt()/BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]", "", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";p_{T}^{after kin. fit}/p_{T}^{gen}",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
    }
    //	  hists.push_back(MyHistogram("B1TrueCorrFactor","BRegJet.genJetPt["+topBranchName_+"recoJetIdxB1]/BRegJet.jetPtCorr["+topBranchName_+"recoJetIdxB1]", "",";p_{T}^{gen}/p_{T}^{corr}; Events",20,0,2 ));
    //	  hists.push_back(MyHistogram("B1GBRFactor","BRegJet.BRegGBRTrainResult["+topBranchName_+"recoJetIdxB1]", "",";regression factor; Events",20,0,2 ));
    //	  hists.back().SetFitGaussToCore();
  }

  if(plotSelectedForPlotting.find("GBRTestingAll")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo;
    helperMyBRegVarInfo.initBRegTest();
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("PlainResponse"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "","; Response (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("PlainResponseB1Only"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0","; Response (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("PlainResponseB1FPGT02Only"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0&&top.fitProb>0.2","; Response (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("PlainResponseB1FPGT02OnlyCP"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0&&top.fitProb>0.2&&top.combinationType==1","; Response (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().SetFitGaussToCore();

      hists.push_back(MyHistogram("FitWRecoB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.fitW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("FitWRecoB1CP"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.fitW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("RecoWRecoB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.recoW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.recoW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.recoW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.recoW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("GenPartonWRecoB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("GenPartonWRecoB1CP"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("GBRFactorB1FPGT02Only"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i),"top.recoB1.Pt()>0&&top.fitProb>0.2","; GBRCorr="+helperMyBRegVarInfo.varNames.at(bvar_i)+"; Events",50,0,2));
      hists.back().SetFitGaussToCore();
    }
    hists.push_back(MyHistogram("GBRFactorUsedForTree","24Cleaned100TreesVarTrainingWithFlatPtTraining","","; regression factor; Events",50,0,2));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("PlainResponseVsMuonFraction","jet_fMuon","24Cleaned100TreesVarTrainingWithFlatPtTraining*BRegJet_jetPtCorr/BRegJet.genJetPt", "",";muon fraction; Response",20,0.01,1,50,0,2 ));
    hists.push_back(MyHistogram("PlainResponseVsMuonFractionBOnly","jet_fMuon","24Cleaned100TreesVarTrainingWithFlatPtTraining*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";muon fraction; Response",20,0.01,1,50,0,2 ));
    hists.push_back(MyHistogram("PlainResponseVsMuonFractionB1Only","jet_fMuon","24Cleaned100TreesVarTrainingWithFlatPtTraining*BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0",";muon fraction; Response",20,0.01,1,50,0,2 ));

    hists.push_back(MyHistogram("PlainResponseDeviationVsJetPtCorrBOnly","BRegJet_jetPtCorr","TMath::Abs(24Cleaned100TreesVarTrainingWithFlatPtTraining-BRegJet.genJetPt/BRegJet_jetPtCorr)", "TMath::Abs(jet.flavour)==5",";jet p_{T}; |GBR-true corr. factor|",20,0.01,100,50,0,2 ));
  }

  if(plotSelectedForPlotting.find("GBRTestingResponseVsInputVar")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo;
    helperMyBRegVarInfo.initBRegInputVarTest();
    std::string BRegUnderStudy("24PlVars1000TreesMin00250EvtsMCS3TQ07");
    std::string BRegUnderStudyName("Nominal");
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("PlainResponseBOnlyVs"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i), BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5", ";"+helperMyBRegVarInfo.varNames.at(bvar_i)+";Response",20,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i), 20, 0, 2));
      hists.back().ConfigurePlotRanges(0.95, 1.05, 0.8, 1.3);

    }
    hists.push_back(MyHistogram("FitWRecoB1"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top}^{fitW} (GBR="+BRegUnderStudyName+"); Events",50,100,250));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("FitWRecoB1CP"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; m_{top}^{fitW} (GBR="+BRegUnderStudyName+"); Events",50,100,250));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("RecoWRecoB1"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"TMath::Sqrt(TMath::Power(top.recoW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.recoW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.recoW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.recoW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top}^{recoW} (GBR="+BRegUnderStudyName+"); Events",50,100,250));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("GenPartonWRecoB1"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+"); Events",50,100,250));
    hists.back().SetFitGaussToCore();
    hists.push_back(MyHistogram("GenPartonWRecoB1CP"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+"); Events",50,100,250));
    hists.back().SetFitGaussToCore();

    hists.push_back(MyHistogram("FitWRecoB1VsJetCorrPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet_jetPtCorr","TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; p_{T}^{jet}; m_{top}^{fitW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));
    hists.push_back(MyHistogram("FitWRecoB1CPVsJetCorrPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet_jetPtCorr","TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; p_{T}^{jet}; m_{top}^{fitW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));

    hists.push_back(MyHistogram("FitWRecoB1VsGenPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet.genJetPt","TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; p_{T}^{gen}; m_{top}^{fitW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));
    hists.push_back(MyHistogram("FitWRecoB1CPVsGenPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet.genJetPt","TMath::Sqrt(TMath::Power(top.fitW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; p_{T}^{gen}; m_{top}^{fitW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));



    hists.push_back(MyHistogram("genpartonWRecoB1VsJetCorrPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet_jetPtCorr","TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; p_{T}^{jet}; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));
    hists.push_back(MyHistogram("genpartonWRecoB1CPVsJetCorrPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet_jetPtCorr","TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; p_{T}^{jet}; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));

    hists.push_back(MyHistogram("genpartonWRecoB1VsGenPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet.genJetPt","TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; p_{T}^{gen}; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));
    hists.push_back(MyHistogram("genpartonWRecoB1CPVsGenPt"+HelperFunctions::cleanedName(BRegUnderStudy.c_str()),"BRegJet.genJetPt","TMath::Sqrt(TMath::Power(top.genpartonW1.E() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.E(),2) - (TMath::Power(top.genpartonW1.Px() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) + "*top.recoB1.Px(),2) + TMath::Power(top.genpartonW1.Py() + " + HelperFunctions::cleanedName(BRegUnderStudy.c_str())+"*top.recoB1.Py(),2) + TMath::Power(top.genpartonW1.Pz() + "+ HelperFunctions::cleanedName(BRegUnderStudy.c_str()) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2&&top.combinationType==1","; p_{T}^{gen}; m_{top}^{genpartonW} (GBR="+BRegUnderStudyName+")",25,30,180,50,100,250));
    helperMyBRegVarInfo.initBRegTest();
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("ResponseCheckBOnlyInclusive"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5","; Response; Events",50,0,2 ));
      hists.back().ConfigureExtraOptions(true);
    }
  }

  if(plotSelectedForPlotting.find("GBRTestingB1Plots")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo;
    helperMyBRegVarInfo.initBRegTest();
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("PlainResponseB1FPGT02Only"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0&&top.fitProb>0.2","; Response (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,helperMyBRegVarInfo.xMins.at(bvar_i), helperMyBRegVarInfo.xMaxs.at(bvar_i) ));
      hists.back().SetFitGaussToCore();

      hists.push_back(MyHistogram("FitWRecoB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.fitW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.fitW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.fitW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.fitW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
      hists.push_back(MyHistogram("RecoWRecoB1"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i)),"TMath::Sqrt(TMath::Power(top.recoW1.E() + " + helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.E(),2) - (TMath::Power(top.recoW1.Px() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) + "*top.recoB1.Px(),2) + TMath::Power(top.recoW1.Py() + " + helperMyBRegVarInfo.varForms.at(bvar_i)+"*top.recoB1.Py(),2) + TMath::Power(top.recoW1.Pz() + "+ helperMyBRegVarInfo.varForms.at(bvar_i) +"*top.recoB1.Pz(),2)))","top.fitProb>0.2","; m_{top} (GBR="+helperMyBRegVarInfo.varNames.at(bvar_i)+"); Events",50,100,250));
      hists.back().SetFitGaussToCore();
    }
  }

  if(plotSelectedForPlotting.find("GBRSystTestingPlots")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo;
    helperMyBRegVarInfo.initBRegTest();
    //	  helperMyBRegVarInfo.initBRegSystematicsTest();
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("PlainResponseBOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      hists.push_back(MyHistogram("PlainResponseBOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      hists.push_back(MyHistogram("PlainResponseNoB"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)!=5",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      hists.push_back(MyHistogram("PlainResponseBTagReq"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      hists.push_back(MyHistogram("PlainResponseBTagReqExtreme"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.9",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
    }
  }


  if(plotSelectedForPlotting.find("GBRTestingResponseBTails")!=plotSelectedForPlotting.end()){
    hists.push_back(MyHistogram("PlainResponseB1TailTest","BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0&&top.fitProb>0.2","; Response; Events",50,0,2));
    hists.push_back(MyHistogram("PlainResponseB1TailTestLogY","BRegJet_jetPtCorr/BRegJet.genJetPt", "top.recoB1.Pt()>0&&top.fitProb>0.2","; Response; Events",50,0,2));
    hists.back().ConfigureExtraOptions(false, "1.",true,"",false,true);

    hists.push_back(MyHistogram("PlainResponseBTailTest","BRegJet_jetPtCorr/BRegJet.genJetPt", "","; Response; Events",50,0,2));
    hists.push_back(MyHistogram("PlainResponseBTailTestLogY","BRegJet_jetPtCorr/BRegJet.genJetPt", "","; Response; Events",50,0,2));
    hists.back().ConfigureExtraOptions(false, "1.",true,"",false,true);
  }

  if(plotSelectedForPlotting.find("CalibTreeResponseBTails")!=plotSelectedForPlotting.end()){
//    hists.push_back(MyHistogram("CalibPlainResponseTailTest","JetPt*JetCorrL1*JetCorrL2L3/GenJetPt", "JetGenJetDeltaR<0.1","; Response; Events",50,0,2));
//    hists.push_back(MyHistogram("CalibPlainResponseTailTestLogY","JetPt*JetCorrL1*JetCorrL2L3/GenJetPt", "JetGenJetDeltaR<0.1","; Response; Events",50,0,2));
//    hists.back().ConfigureExtraOptions(false, "1.",true,"",false,true);

    hists.push_back(MyHistogram("CalibPlainResponseLeadJetTailTest","JetPt[0]*JetCorrL1[0]*JetCorrL2L3[0]/GenJetPt[0]", "JetGenJetDeltaR[0]<0.1","; Response; Events",50,0,2));
    hists.push_back(MyHistogram("CalibPlainResponseLeadJetTailTestLogY","JetPt[0]*JetCorrL1[0]*JetCorrL2L3[0]/GenJetPt[0]", "JetGenJetDeltaR[0]<0.1","; Response; Events",50,0,2));
    hists.back().ConfigureExtraOptions(false, "1.",true,"",false,true);
  }


  if(plotSelectedForPlotting.find("GBRResponseTesting")!=plotSelectedForPlotting.end()){
    MyBRegVarInfo helperMyBRegVarInfo;
    helperMyBRegVarInfo.initBRegSystematicsTest();
    //	  helperMyBRegVarInfo.initBRegTest();
    for (size_t bvar_i=0;bvar_i<helperMyBRegVarInfo.varNames.size();++bvar_i){
      hists.push_back(MyHistogram("PlainResponseBOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      //	    void TopMassControlPlots::MyHistogram::ConfigureExtraOptions(bool SetFitGaussToCore, std::string SetCustomHistoweight, bool ExportSigVarToRoot, std::string LegendHeader, bool useLogXBins, bool useLogYBins)

      hists.push_back(MyHistogram("PlainResponseBOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);

      hists.push_back(MyHistogram("PlainResponseUDSOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)<=3",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);
      hists.push_back(MyHistogram("PlainResponseUDSOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)<=3",";jet p_{T}^{gen}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);


      hists.push_back(MyHistogram("PlainResponseCOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==4",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);

      hists.push_back(MyHistogram("PlainResponseNoB"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)!=5",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);

      hists.push_back(MyHistogram("PlainResponseBTagReq"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);

      hists.push_back(MyHistogram("PlainResponseBTagReqExtreme"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsJetPtCorr","BRegJet_jetPtCorr",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.9",";jet p_{T}; Response",20,30.,300,50,0,2 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);


      hists.push_back(MyHistogram("RecoPtLogLog"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",20,30.,300,50,30,300));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true,true);

      hists.push_back(MyHistogram("RecoPtLog"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",20,30.,300,50,30,300));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true,false);

      hists.push_back(MyHistogram("FineBinBOnlyRecoPtLog"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",200,5,300,100,30,500));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true,false);
      hists.push_back(MyHistogram("FineBinPlainResponseBOnly"+HelperFunctions::cleanedName(helperMyBRegVarInfo.varNames.at(bvar_i))+"VsGenJetPt","BRegJet.genJetPt",helperMyBRegVarInfo.varForms.at(bvar_i)+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5",";p_{T}^{gen}; Response",200,5,300,100,0,4 ));
      hists.back().ConfigureExtraOptions(false, "1.", true, "", true);



    }
  }

  if(plotSelectedForPlotting.find("GBRResponseImprovement")!=plotSelectedForPlotting.end()){
//	  std::string BRegUnderStudy("24PlVars1000TreesMin00250EvtsMCS3TQ07");

	  std::string BRegUnderStudy("24Cleaned100TreesVarTrainingWithFlatPtTraining");
	  std::cout << "adding GBRResponseImprovement histos" << std::endl;
	  hists.push_back(MyHistogram("JetPtCorr",BRegUnderStudy+"*BRegJet_jetPtCorr", "","; p_{T}^{corr.}; Events",20,20,500 ));
	  hists.back().ConfigureExtraOptions(false, "1.", false,"",true);
//	  hists.back().ConfigureExtraOptions(false, "1.", false,"",false);
//
	  hists.push_back(MyHistogram("ResponseCheckBOnlyInclusive",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"");
//
//
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyLowPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>30&&BRegJet_jetPtCorr<50","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; p_{T}<50 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyPeakPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>50&&BRegJet_jetPtCorr<70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; 50<p_{T}<70 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyHighPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; p_{T}>70 GeV");
//
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyLowPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>30&&BRegJet_jetPtCorr<50","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; p_{T}<50 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyPeakPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>50&&BRegJet_jetPtCorr<70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; 50<p_{T}<70 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBOnlyHighPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "TMath::Abs(jet.flavour)==5&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; p_{T}>70 GeV");
//
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqLowPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>30&&BRegJet_jetPtCorr<50","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; p_{T}<50 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqPeakPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>50&&BRegJet_jetPtCorr<70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; 50<p_{T}<70 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqHighPtBarrel",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)<1.3&&BRegJet_jetPtCorr>70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|<1.3; p_{T}>70 GeV");
//
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqLowPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>30&&BRegJet_jetPtCorr<50","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; p_{T}<50 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqPeakPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>50&&BRegJet_jetPtCorr<70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; 50<p_{T}<70 GeV");
//	  hists.push_back(MyHistogram("ResponseCheckBTagReqHighPtEndcap",BRegUnderStudy+"*BRegJet_jetPtCorr/BRegJet.genJetPt", "jet_bTagCSV>0.679&&TMath::Abs(BRegJet_jetEta)>1.3&&BRegJet_jetPtCorr>70","; Response; Events",50,0,2 ));
//	  hists.back().ConfigureExtraOptions(true, "1.", false,"|#eta|>1.3; p_{T}>70 GeV");
//
//
//
//	  hists.push_back(MyHistogram("PlainResponseDeviationVsJetPtCorrBOnly","BRegJet_jetPtCorr","TMath::Abs("+BRegUnderStudy+"-BRegJet.genJetPt/BRegJet_jetPtCorr)", "TMath::Abs(jet.flavour)==5",";jet p_{T}; |GBR-true corr. factor|",20,30,300,50,0,2 ));
//	  hists.back().ConfigureExtraOptions(false, "1.", false,"",true);
//	  hists.back().ConfigurePlotRanges(0.7, 1.0, 0.0, 0.3);
//	  hists.push_back(MyHistogram("PlainResponseDeviationVsJetEtaBOnly","TMath::Abs(BRegJet_jetEta)","TMath::Abs("+BRegUnderStudy+"-BRegJet.genJetPt/BRegJet_jetPtCorr)", "TMath::Abs(jet.flavour)==5",";jet |#eta|; |GBR-true corr. factor|",12,0,2.4,50,0,2 ));
//	  hists.back().ConfigurePlotRanges(0.8, 1.0,0.0,0.25);

	  std::cout << "done adding GBRResponseImprovement histos" << std::endl;
  }
  //
  // DEFINE DATASETS
  //
  
  std::cout << "define datasets" << std::endl;
  // Alljets channel
  if (channelID == Helper::kAllJets) {
    samples.push_back(MySample("Data", "Run2012_alljets*", kData, kBlack));
    samples.push_back(MySample("t#bar{t}", "Summer12_TTJetsMS1725_1.00_alljets*", kSig, kRed+1, 1, 1., "", 0.04/po::GetOption<double>("templates.fSig")));

    if(plotSelectedForPlotting.find("MassVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("m_{t} = 169.5 GeV", "Summer12_TTJetsMS1695_1.00_alljets*", kSigVar, kBlue+1 , 3));
      samples.push_back(MySample("m_{t} = 175.5 GeV", "Summer12_TTJetsMS1755_1.00_alljets*", kSigVar, kGreen+1, 2)); 
    }

    // JES
    if(plotSelectedForPlotting.find("JESUncVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, JES down", "Summer12_TTJetsMS1725_source:down_Total_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t}, JES up", "Summer12_TTJetsMS1725_source:up_Total_alljets*", kSigVar, kGreen+1, 1));
    }
    
    // JER
    if(plotSelectedForPlotting.find("JERUncVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, JER down", "Summer12_TTJetsMS1725_jer:-1_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t}, JER up", "Summer12_TTJetsMS1725_jer:1_alljets*", kSigVar, kGreen+1, 1));
    }    
    
    // PU
    if(plotSelectedForPlotting.find("PUPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, PU up", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kRed+1, 1, 1., "weight.combinedWeight|weight.combinedWeight/weight.puWeight*weight.puWeightUp"));
      samples.push_back(MySample("t#bar{t}, PU down", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kRed+1, 1, 1., "weight.combinedWeight|weight.combinedWeight/weight.puWeight*weight.puWeightDown"));
    }
    
    // MCGeneratorPlots
    if(plotSelectedForPlotting.find("MCGeneratorPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg_alljets*", kSigVar, kGreen+1, 9));
      samples.push_back(MySample("t#bar{t} scaleup", "Summer12_TTJetsMS1725_scaleup_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t} scaledown", "Summer12_TTJetsMS1725_scaledown_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t} matchingup", "Summer12_TTJetsMS1725_matchingup_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t} matchingdown", "Summer12_TTJetsMS1725_matchingdown_alljets*", kSigVar, kRed+1, 1));
      samples.push_back(MySample("t#bar{t}, top pt reweighting", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kBlack, 1, 1., "weight.combinedWeight|weight.combinedWeight*sqrt(exp(0.156-0.00137*top.genpartonTop1.Pt())*exp(0.156-0.00137*top.genpartonTop2.Pt()))"));
    }

    // TriggerPlots
    if(plotSelectedForPlotting.find("TriggerPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, triggerup", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kBlack, 1, 1.  , "(0.5*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.5)|(0.25*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.75)"));
      samples.push_back(MySample("t#bar{t}, triggerdown", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kBlack, 1, 1., "(0.5*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.5)|(0.75*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.25)"));
    }

    // bJES
    if(plotSelectedForPlotting.find("BJESPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, frag", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kRed+1, 1, 1., "weight.combinedWeight|weight.combinedWeight*(weight.bJESWeight_frag==-1?1:weight.bJESWeight_frag)"));
      samples.push_back(MySample("t#bar{t}, #nu up", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kRed+1, 1, 1., "weight.combinedWeight|weight.combinedWeight*(weight.bJESWeight_frag==-1?1:weight.bJESWeight_fNuUp)"));
      samples.push_back(MySample("t#bar{t}, #nu down", "Summer12_TTJetsMS1725_1.00_alljets*", kSigVar, kRed+1, 1, 1., "weight.combinedWeight|weight.combinedWeight*(weight.bJESWeight_frag==-1?1:weight.bJESWeight_fNuDown)"));
    }
    
    /*
    // Signal variations (UE Tune)
    if(plotSelectedForPlotting.find("UETunePlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, Z2*"     , "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kRed+1    , 1));
      samples.push_back(MySample("t#bar{t}, P11"     , "Z2_S12_P11_sig"              , kSigVar, kMagenta+1, 9));
      samples.push_back(MySample("t#bar{t}, P11mpiHi", "Z2_S12_P11mpiHi_sig"         , kSigVar, kBlue+1   , 3));
      samples.push_back(MySample("t#bar{t}, P11TeV"  , "Z2_S12_P11TeV_sig"           , kSigVar, kGreen+4  , 4));
      samples.push_back(MySample("t#bar{t}, P11noCR" , "Z2_S12_P11NoCR_sig"          , kSigVar, kCyan+1   , 2));
    }

    // Signal variations (MC Generator)
    if(plotSelectedForPlotting.find("MCGeneratorPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, MG Z2* MadSpin", "Z2_S12_ABS_JES_100_172_5_MadSpin_sig*", kSigVar, kRed+1  , 1));
      samples.push_back(MySample("t#bar{t}, MG Z2*"        , "Z2_S12_ABS_JES_100_172_5_sig", kSigVar, kRed  , 1));
      samples.push_back(MySample("t#bar{t}, Powheg+Pythia", "Z2_S12_POWHEG_sig"           , kSigVar, kGreen+1, 9));
      samples.push_back(MySample("t#bar{t}, Powheg+Herwig", "Z2_S12_POWHER_sig"           , kSigVar, kBlue+1 , 7));
      samples.push_back(MySample("t#bar{t}, MC@NLO"       , "Z2_S12_MCNLO_sig"            , kSigVar, kCyan+1 , 2));
    }

    // Signal variations (JES)
    if(plotSelectedForPlotting.find("JESVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, JES = 0.96", "Z2_S12_ABS_JES_096_172_5_MadSpin_sig", kSigVar, kRed+1    , 1));
      samples.push_back(MySample("t#bar{t}, JES = 0.98", "Z2_S12_ABS_JES_098_172_5_MadSpin_sig", kSigVar, kMagenta+1, 1));
      samples.push_back(MySample("t#bar{t}, JES = 1.00", "Z2_S12_ABS_JES_100_172_5_MadSpin_sig", kSigVar, kBlue+1   , 1));
      samples.push_back(MySample("t#bar{t}, JES = 1.02", "Z2_S12_ABS_JES_102_172_5_MadSpin_sig", kSigVar, kCyan+1   , 1));
      samples.push_back(MySample("t#bar{t}, JES = 1.04", "Z2_S12_ABS_JES_104_172_5_MadSpin_sig", kSigVar, kGreen+1  , 1));
    }
    */
    samples.push_back(MySample("Background", "Background_MJP12", kBkg, kYellow, 1));
    //samples.push_back(MySample("Background", "Run2012_Mixing8_alljets*", kBkg, kYellow, 1)); 
    //samples.push_back(MySample("Background", "Run2012_Mixing8_fix_alljets*", kBkg, kYellow, 1)); 
    //samples.push_back(MySample("Background", "job_QCDMixing_MJPS12*", kBkg, kYellow, 1)); 
    samples.push_back(MySample("Background (mixing)", "mix6_QCDMixing_MJPS12*", kBkgVar, kYellow, 1)); 
  }
  

  else if(plotSelectedForPlotting.find("CalibTreeMaker")!=plotSelectedForPlotting.end()){
	  std::cout << "going to add CalibTreeMaker samples" << std::endl;
	  if(plotSelectedForPlotting.find("PYTHIAQCD")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("", "PYTHIA6", kData, kWhite));
	    samples.push_back(MySample("PYTHIA QCD", "PYTHIA6", kSig, kRed+1));
	  }

  }

  else if(plotSelectedForPlotting.find("GBRTesting")!=plotSelectedForPlotting.end()){
	  std::cout << "going to add GBRTesting samples" << std::endl;
	  if(plotSelectedForPlotting.find("BasicTests")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1/0.0359779637));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1, 1, 1/0.0359779637));
//	    samples.push_back(MySample("t#bar{t}(1.04)", "MC104NoBReg", kSigVar, kAzure-2));
////	    samples.push_back(MySample("t#bar{t}(1.04, BReg)", "MC104", kSigVar, kGreen-3));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSigVar, kCyan+1));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO, BReg)", "MCatNLO", kSigVar, kMagenta));
//	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
//	    samples.push_back(MySample("t#bar{t}(POWHER, BReg)", "POWHER", kSigVar, kGreen-5));
	    samples.push_back(MySample("Data", "DataNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("Data (BReg)", "Data", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1, 1, 1/0.0359779637));
	  }
	  else if(plotSelectedForPlotting.find("BasicTestsNoData")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
//	    samples.push_back(MySample("t#bar{t}(1.04)", "MC104NoBReg", kSigVar, kAzure-2));
////	    samples.push_back(MySample("t#bar{t}(1.04, BReg)", "MC104", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(MC@NLO, BReg)", "MCatNLO", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(POWHER, BReg)", "POWHER", kSigVar, kGreen-5));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("BasicTestsNoBRegNoData")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kData, kBlack,1,1/0.0359779637));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSig, kRed));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSigVar, kRed));
//	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSig, kRed));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kRed));
	    samples.push_back(MySample("t#bar{t}, #nu up", "MCNoBReg", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
	    samples.push_back(MySample("t#bar{t}, frag up", "MCNoBReg", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
	  }
	  else if(plotSelectedForPlotting.find("BasicTestsBRegNoData")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t}", "MC", kData, kBlack,1,1/0.0359779637));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLO", kSig, kRed));
//	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLO", kSigVar, kRed));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHER", kSig, kRed));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHER", kSigVar, kRed));
	    samples.push_back(MySample("t#bar{t}, #nu up", "MC", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
	    samples.push_back(MySample("t#bar{t}, frag", "MC", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
	  }

	  else if(plotSelectedForPlotting.find("BasicTestsAllVariationsNoData")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kData, kBlack,1,1));
	    samples.push_back(MySample("t#bar{t}(MadSpin)", "MadSpinNoBReg", kSig, kRed)); 

	    samples.push_back(MySample("t#bar{t}(MadSpin)", "MadSpinNoBReg", kSigVar, kRed));
	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSigVar, kBlack)); 
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(POWHEG)", "POWHEGNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}, #nu up", "MadSpinNoBReg", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
	    samples.push_back(MySample("t#bar{t}, #nu down", "MadSpinNoBReg", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuDown"));
	    samples.push_back(MySample("t#bar{t}, frag", "MadSpinNoBReg", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
	    samples.push_back(MySample("t#bar{t}, fragHard", "MadSpinNoBReg", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fragHard"));
	    samples.push_back(MySample("t#bar{t}, fragSoft", "MadSpinNoBReg", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fragSoft"));
	    samples.push_back(MySample("t#bar{t}(BUp)", "FlavPureBUpNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(BDn)", "FlavPureBDownNoBReg", kSigVar, kCyan-4));


	    samples.push_back(MySample("t#bar{t}BReg(MadSpin)", "MadSpin", kSigVar, kRed));
	    samples.push_back(MySample("t#bar{t}BReg(MC@NLO)", "MCatNLO", kSigVar, kBlack)); 
	    samples.push_back(MySample("t#bar{t}BReg(POWHER)", "POWHER", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}BReg(POWHEG)", "POWHEG", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}BReg, #nu up", "MadSpin", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
	    samples.push_back(MySample("t#bar{t}BReg, #nu down", "MadSpin", kSigVar, kMagenta+1, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuDown"));
	    samples.push_back(MySample("t#bar{t}BReg, frag", "MadSpin", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
	    samples.push_back(MySample("t#bar{t}BReg, fragHard", "MadSpin", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fragHard"));
	    samples.push_back(MySample("t#bar{t}BReg, fragSoft", "MadSpin", kSigVar, kGreen-5, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fragSoft"));
	    samples.push_back(MySample("t#bar{t}BReg(BUp)", "FlavPureBUp", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}BReg(BDn)", "FlavPureBDown", kSigVar, kCyan-4));

	  }

	  else if(plotSelectedForPlotting.find("FragmentationTestsNoData")!=plotSelectedForPlotting.end()){
//		    samples.push_back(MySample("t#bar{t}BReg", "MadSpin", kData, kBlack));
//		    samples.push_back(MySample("t#bar{t}BReg, frag", "MadSpin", kSig, kRed, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
//		    samples.push_back(MySample("t#bar{t}BReg, frag", "MadSpin", kSigVar, kRed, 1, 1,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));

		    samples.push_back(MySample("t#bar{t}BReg", "MC", kData, kBlack, 1, 1/0.0359779637));
		    samples.push_back(MySample("t#bar{t}BReg, frag", "MC", kSig, kRed, 1, 1/0.0359779637,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
		    samples.push_back(MySample("t#bar{t}BReg, frag", "MC", kSigVar, kRed, 1, 1/0.0359779637,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));

	  }

	  else if(plotSelectedForPlotting.find("ExtremeBasicTestsNoData")!=plotSelectedForPlotting.end()){
		  std::cout << "adding ExtremeBasicTestsNoData samples" << std::endl;
//		    MySample(std::string name_, std::string file_, int type_, int color_, int line_ = 1, double scale_ = 1., std::string replaceVar_="") :

		    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kData, kBlack, 1, 1/0.0359779637)); //scale with 1/MCweight
	//	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
		    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kSig, kRed+1, 1, 1/0.0359779637));

		    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kSigVar, kRed+1, 1, 1/0.0359779637));

//		    samples.push_back(MySample("t#bar{t}", "MadSpinNoBReg", kData, kBlack, 1, 1/0.0040092475));
//		    samples.push_back(MySample("t#bar{t} (BReg)", "MadSpin", kSig, kRed+1, 1, 1/0.0040092475));
//
//		    samples.push_back(MySample("t#bar{t} (BReg)", "MadSpin", kSigVar, kRed+1, 1, 1/0.0040092475));
	  }
	  else if(plotSelectedForPlotting.find("ExtremeBasicTestsNoDataHideBReg")!=plotSelectedForPlotting.end()){
		  std::cout << "adding ExtremeBasicTestsNoData samples" << std::endl;
		    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kData, kBlack, 1, 1/0.0359779637)); //scale with 1/MCweight
//		    samples.push_back(MySample("", "MC", kSig, kWhite, 1, 1/0.0359779637));
//		    samples.push_back(MySample("", "MC", kSigVar, kWhite, 1, 1/0.0359779637));

	  }
	  else if(plotSelectedForPlotting.find("VeryBasicTestsNoData")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));

	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(POWHER, BReg)", "POWHER", kSigVar, kGreen-5));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("JESUncTests")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(BUp)", "FlavPureBUpNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(BUp, BReg)", "FlavPureBUp", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(BDn)", "FlavPureBDownNoBReg", kSigVar, kCyan-4));
	    samples.push_back(MySample("t#bar{t}(BDn, BReg)", "FlavPureBDown", kSigVar, kMagenta-4));
	    samples.push_back(MySample("t#bar{t}(JESUp)", "JESUncUpNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(JESUp, BReg)", "JESUncUp", kSigVar, kGreen-5));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("GlobalJES")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(1.04)", "MC104NoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(1.04, BReg)", "MC104", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(0.96)", "MC096NoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(0.96, BReg)", "MC096", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("MassDiff")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(178.5)", "MC1785NoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(175.5, BReg)", "MC1785", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(166.5)", "MC1665NoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(166.5, BReg)", "MC1665", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("MatchingTests")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(MUp)", "MatchingupNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(MUp, BReg)", "Matchingup", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(MDn)", "MatchingDownNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(MDn, BReg)", "MatchingDown", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("ScaleTests")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(SUp)", "ScaleupNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(SUp, BReg)", "Scaleup", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(SDn)", "ScaleDownNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(SDn, BReg)", "ScaleDown", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
	  else if(plotSelectedForPlotting.find("AllGBRTests")!=plotSelectedForPlotting.end()){
	    samples.push_back(MySample("t#bar{t} (BReg)", "MC", kData, kBlack, 1, 1));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSig, kRed+1));
	    samples.push_back(MySample("t#bar{t}(MC@NLO)", "MCatNLONoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(MC@NLO, BReg)", "MCatNLO", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(POWHER)", "POWHERNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(POWHER, BReg)", "POWHER", kSigVar, kGreen-5));
	    samples.push_back(MySample("Data", "DataNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("Data, BReg)", "Data", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(BUp)", "FlavPureBUpNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(BUp, BReg)", "FlavPureBUp", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(BDn)", "FlavPureBDownNoBReg", kSigVar, kCyan+5));
	    samples.push_back(MySample("t#bar{t}(BDn, BReg)", "FlavPureBDown", kSigVar, kMagenta+5));
	    samples.push_back(MySample("t#bar{t}(JESUp)", "JESUncUpNoBReg", kSigVar, kAzure-5));
	    samples.push_back(MySample("t#bar{t}(JESUp, BReg)", "JESUncUp", kSigVar, kGreen-5));
	    samples.push_back(MySample("t#bar{t}(1.04)", "MC104NoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(1.04, BReg)", "MC104", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(0.96)", "MC096NoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(0.96, BReg)", "MC096", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(178.5)", "MC1785NoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(175.5, BReg)", "MC1785", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(166.5)", "MC1665NoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(166.5, BReg)", "MC1665", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(MUp)", "MatchingupNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(MUp, BReg)", "Matchingup", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(MDn)", "MatchingDownNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(MDn, BReg)", "MatchingDown", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}(SUp)", "ScaleupNoBReg", kSigVar, kAzure-2));
	    samples.push_back(MySample("t#bar{t}(SUp, BReg)", "Scaleup", kSigVar, kGreen-3));
	    samples.push_back(MySample("t#bar{t}(SDn)", "ScaleDownNoBReg", kSigVar, kCyan+1));
	    samples.push_back(MySample("t#bar{t}(SDn, BReg)", "ScaleDown", kSigVar, kMagenta));
	    samples.push_back(MySample("t#bar{t}", "MCNoBReg", kSigVar, kRed+1));
	  }
  }

  // Lepton+jets channel
  else {
    // DATA
    if(plotSelectedForPlotting.find("ThesisData")!=plotSelectedForPlotting.end()) {
      samples.push_back(MySample("Data", "Run2015D", kData, kBlack));
    }
    else {
        if (channelID != Helper::kElectronJets) samples.push_back(MySample("Data", "76x/SingleMuon_16Dec2015_workingJEC_allCut", kData, kBlack));
    	//if (channelID != Helper::kElectronJets) samples.push_back(MySample("Data", "76x/RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_workingJEC_1.00_allCut", kData, kBlack,1,lumi_/1000.));

    	if (channelID != Helper::kMuonJets) samples.push_back(MySample("Data", "76x/SingleElectron_16Dec2015_workingJEC_allCut", kData, kBlack));
    }
    
    // SIGNAL
    double signalNormUnc = 0.0;
    if(plotSelectedForPlotting.find("BackgroundPlots")!=plotSelectedForPlotting.end()) signalNormUnc = 0.053; //TODO?
    
    if(plotSelectedForPlotting.find("NobPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}", "Summer12_TTJetsMS1725_nob", kSig, kRed+1, 1, lumi_/1000., "", signalNormUnc));
    }
    else samples.push_back(MySample("t#bar{t}", "76x/RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_sysJER_1.00_sigRes0.0_allCut", kSig, kRed+1, 1, lumi_/1000., "", signalNormUnc));
    	//samples.push_back(MySample("t#bar{t}", "76x/RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_workingJEC_1.00_allCut", kSig, kRed+1, 1, lumi_/1000., "", signalNormUnc));
    //samples.push_back(MySample("t#bar{t}", "../2015_JESVariations/76x/RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_patJetsUpdated_1.00", kSig, kRed+1, 1, lumi_/1000., "", signalNormUnc));

    // SIGNAL VARIATIONS
    
    // JES
    if(plotSelectedForPlotting.find("JESUncVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, JES down", "Summer12_TTJetsMS1725_source:down_Total", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JES up", "Summer12_TTJetsMS1725_source:up_Total", kSigVar, kGreen+1, 1, lumi_/1000.));
    }
    
    // JER
    if(plotSelectedForPlotting.find("JERUncVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, JER down", "Summer12_TTJetsMS1725_jer:-1", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JER up", "Summer12_TTJetsMS1725_jer:1", kSigVar, kGreen+1, 1, lumi_/1000.));
    }    
    
    // PU
    if(plotSelectedForPlotting.find("PUPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, PU up", "Summer12_TTJetsMS1725_nob", kSigVar, kRed+1, 1, lumi_/1000., "weight.combinedWeight|weight.combinedWeight/weight.puWeight*weight.puWeightUp"));
      samples.push_back(MySample("t#bar{t}, PU down", "Summer12_TTJetsMS1725_nob", kSigVar, kRed+1, 1, lumi_/1000., "weight.combinedWeight|weight.combinedWeight/weight.puWeight*weight.puWeightDown"));
    }
    
    // MCGeneratorPlots
    if(plotSelectedForPlotting.find("MCGeneratorPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} scaleup", "Summer12_TTJetsMS1725_scaleup", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} scaledown", "Summer12_TTJetsMS1725_scaledown", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} matchingup", "Summer12_TTJetsMS1725_matchingup", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} matchingdown", "Summer12_TTJetsMS1725_matchingdown", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, top pt reweighting", "Summer12_TTJetsMS1725_1.00", kSigVar, kBlack, 1, lumi_/1000., "weight.combinedWeight|weight.combinedWeight*sqrt(exp(0.156-0.00137*top.genpartonTop1.Pt())*exp(0.156-0.00137*top.genpartonTop2.Pt()))"));
    }
    
    // MassPlots
    if(plotSelectedForPlotting.find("MassPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t} 1715", "Summer12_TTJetsMS1715_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} 1735", "Summer12_TTJetsMS1735_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
    }
    
    // ShowerHadronizationPlots
    if(plotSelectedForPlotting.find("ShowerHadronizationPlots")!=plotSelectedForPlotting.end()){
      //samples.push_back(MySample("t#bar{t}, MadGraph+Pythia Z2*", "Summer12_TTJetsMS1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} scale up", "Summer12_TTJetsMS1725_scaleup", kSigVar, kRed-1, 5, lumi_/1000.));
      samples.push_back(MySample("t#bar{t} scale down", "Summer12_TTJetsMS1725_scaledown", kSigVar,kRed-7, 6, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Powheg+Herwig6", "Summer12_TTJets1725_powheg_herwig", kSigVar, kOrange+2, 7, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, MC@NLO+Herwig6", "Summer12_TTJets1725_mcatnlo_herwig", kSigVar, kBlue+1, 2, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Sherpa", "Summer12_TTJets1725_sherpa", kSigVar, kYellow+1, 3, lumi_/1000./0.7719));
    }
    
    // ShowerHadronizationFastPlots
    if(plotSelectedForPlotting.find("ShowerHadronizationFastPlots")!=plotSelectedForPlotting.end()){
      //samples.push_back(MySample("t#bar{t}, MadGraph+Pythia Z2*", "Summer12_TTJetsMS1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      //samples.push_back(MySample("t#bar{t} scale up", "Summer12_TTJetsMS1725_scaleup", kSigVar, kRed-1, 5, lumi_/1000.));
      //samples.push_back(MySample("t#bar{t} scale down", "Summer12_TTJetsMS1725_scaledown", kSigVar,kRed-7, 6, lumi_/1000.));
      samples.push_back(MySample("FastSim, t#bar{t}, MadGraph+Pythia Z2*", "Summer12_TTJets1725_FSIM_1.00", kSigVar, kBlack, 1, lumi_/1000.));
      samples.push_back(MySample("FastSim, t#bar{t}, aMC@NLO+Pythia8", "Summer12_TT1725_amcatnlo_pythia8", kSigVar, kGreen+1, 9, lumi_/1000.));
      samples.push_back(MySample("FastSim, t#bar{t}, aMC@NLO+Herwig++", "Summer12_TT1725_amcatnlo_herwigpp", kSigVar, kBlue+1, 2, lumi_/1000.));
      samples.push_back(MySample("FastSim, t#bar{t}, Sherpa Cluster", "Summer12_TT1725_sherpa2_cluster", kSigVar, kYellow+1, 3, lumi_/1000.));
      samples.push_back(MySample("FastSim, t#bar{t}, Sherpa Lund", "Summer12_TT1725_sherpa2_lund", kSigVar, kYellow-5, 3, lumi_/1000.));
    }
    
    // PowhegPlots
    if(plotSelectedForPlotting.find("PowhegPlots")!=plotSelectedForPlotting.end()){
      //samples.push_back(MySample("t#bar{t}, MadGraph+Pythia Z2*", "Summer12_TTJetsMS1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      //samples.push_back(MySample("t#bar{t} scale up", "Summer12_TTJetsMS1725_scaleup", kSigVar, kRed-1, 5, lumi_/1000.));
      //samples.push_back(MySample("t#bar{t} scale down", "Summer12_TTJetsMS1725_scaledown", kSigVar,kRed-7, 6, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Powheg2+Pythia6", "Summer12_TTJets1725_powheg2_pythia6_merge", kSigVar, kOrange+2, 7, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Powheg2+Pythia6 P11C", "Summer12_TTJets1725_powheg_P11C", kSigVar, kBlue+2, 1, lumi_/1000.));
    }
    
    // BACKGROUND
    if(plotSelectedForPlotting.find("BackgroundPlots")!=plotSelectedForPlotting.end()){
	samples.push_back(MySample("Single t", "76x/BG/major/ST_merged_allCut", kBkg, kMagenta, 1, lumi_/1000., "", 0.1));
      samples.push_back(MySample("W+jets", "76x/BG/major/WJetsToLNu_TuneCUETP8M-amcatnloFXFX-pythia8_asymptotic_v12-v1_allCut", kBkg, kGreen-3, 1, lumi_/1000.  , "", 0.2));   
      samples.push_back(MySample("Z+jets", "76x/BG/minor/DYJetsToLL_merged_allCut", kBkg, kAzure-2, 1, lumi_/1000., "", 0.2)); //FIXME check if neg. Weight of amcatnloFXFX part were put into acount
      samples.push_back(MySample("QCD multijet", "76x/BG/minor/QCD_lepEnriched_merged_allCut", kBkg, kYellow, 1, lumi_/1000., "", 1.));
      samples.push_back(MySample("Diboson", "76x/BG/minor/DiBoson_merged_allCut", kBkg, kWhite, 1, lumi_/1000., "", 0.05));
    }

    // OTHER VARIATIONS
    if(plotSelectedForPlotting.find("BasicTestsNoData")!=plotSelectedForPlotting.end()){
    	samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kData, kBlack, 1, 1));
        samples.push_back(MySample("t#bar{t}, MG+Pythia6 P11", "Summer12_TTJets1725_MGDecays_P11", kSig, kMagenta+1, 9, lumi_/1000.*9./4.));
        samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, Powheg+Herwig6 AUET2", "Summer12_TTJets1725_powheg_herwig", kSigVar, kBlue+1, 7, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, MC@NLO", "Summer12_TTJets1725_mcatnlo_herwig", kSigVar, kCyan+1, 1, lumi_/1000.));
    }

    if(plotSelectedForPlotting.find("BasicTests")!=plotSelectedForPlotting.end()){
        samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, Powheg+Pythia6 Z2*", "Summer12_TTJets1725_powheg", kSigVar, kGreen+1, 9, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, Powheg+Herwig6 AUET2", "Summer12_TTJets1725_powheg_herwig", kSigVar, kBlue+1, 7, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, MC@NLO", "Summer12_TTJets1725_mcatnlo_herwig", kSigVar, kCyan+1, 1, lumi_/1000.));
    }

    // UETunePlots
    if(plotSelectedForPlotting.find("UETunePlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, MG+Pythia6 Z2*", "Summer12_TTJets1725_MGDecays", kSigVar, kRed+1, 1, lumi_/1000.*9./4.));
      samples.push_back(MySample("t#bar{t}, MG+Pythia6 P11", "Summer12_TTJets1725_MGDecays_P11", kSigVar, kMagenta+1, 9, lumi_/1000.*9./4.));
      samples.push_back(MySample("t#bar{t}, MG+Pythia6 P11noCR", "Summer12_TTJets1725_MGDecays_P11noCR", kSigVar, kCyan+1, 2, lumi_/1000.*9./4.));
    }

    // SherpaPlots
    if(plotSelectedForPlotting.find("SherpaPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, MG+Pythia6 Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JES = 1.02", "Summer12_TTJets1725_1.02", kSigVar, kMagenta+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Sherpa (4jets)", "Summer12_TTJets1725_sherpa", kSigVar,  kCyan+1, 2, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Sherpa Cluster", "Summer12_TT_cluster_8TeV-sherpa", kSigVar, kGreen+1, 9, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Sherpa Lund", "Summer12_TT_lund_8TeV-sherpa", kSigVar, kBlue+1, 7, lumi_/1000.));
    }

    // JES
    if(plotSelectedForPlotting.find("JESVariationPlots")!=plotSelectedForPlotting.end()){
      //samples.push_back(MySample("t#bar{t}, JES = 0.96", "Summer12_TTJets1725_0.96", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JES = 0.98", "Summer12_TTJets1725_0.98", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JES = 1.00", "Summer12_TTJets1725_1.00", kSigVar, kBlue+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, JES = 1.02", "Summer12_TTJets1725_1.02", kSigVar, kGreen+1, 1, lumi_/1000.));
      //samples.push_back(MySample("t#bar{t}, JES = 1.04", "Summer12_TTJets1725_1.04", kSigVar, kGreen+1, 1, lumi_/1000.));
    }
    
    
    
    if(plotSelectedForPlotting.find("bJESUncVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, bJES down (FlavorQCD)", "Summer12_TTJets1725_flavor:down_FlavorQCD", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, bJES default", "Summer12_TTJets1725_1.00", kSigVar, kBlue+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, bJES up (FlavorQCD)", "Summer12_TTJets1725_flavor:up_FlavorQCD", kSigVar, kGreen+1, 1, lumi_/1000.));
    }

    // Nu fraction and frag function weighting
    if(plotSelectedForPlotting.find("WeightVariationPlots")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kData, kBlack, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, frag", "Summer12_TTJets1725_1.00", kSig, kRed, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
//        MySample(std::string name_, std::string file_, int type_, int color_, int line_ = 1, double scale_ = 1., std::string replaceVar_="") :
//      samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, #nu up", "Summer12_TTJets1725_1.00", kSigVar, kMagenta+1, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
      samples.push_back(MySample("t#bar{t}, #nu down", "Summer12_TTJets1725_1.00", kSigVar, kBlue+1, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuDown"));
      samples.push_back(MySample("t#bar{t}, fragHard", "Summer12_TTJets1725_1.00", kSigVar, kCyan+1, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fragHard"));
    }

    // Nu fraction and frag function weighting
    if(plotSelectedForPlotting.find("BRegSysPlots")!=plotSelectedForPlotting.end()){
        samples.push_back(MySample("t#bar{t}, Powheg+Herwig6 AUET2", "Summer12_TTJets1725_powheg_herwig", kSigVar, kRed, 1, lumi_/1000.));
        samples.push_back(MySample("t#bar{t}, #nu up", "Summer12_TTJets1725_1.00", kSigVar, kMagenta+1, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_fNuUp"));
        samples.push_back(MySample("t#bar{t}, frag", "Summer12_TTJets1725_1.00", kSigVar, kGreen-5, 1, lumi_/1000.,"weight.combinedWeight|weight.combinedWeight*weight.bJESWeight_frag"));
        samples.push_back(MySample("t#bar{t}", "Summer12_TTJets1725_1.00", kSigVar, kBlue+1, 1, lumi_/1000.));
    }
    //Before/after b regression comparison
    if(plotSelectedForPlotting.find("BeforeAfter")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Z2*, noBReg", "Summer12_TTJets1725_1.00", kSigVar, kCyan+1, 1, lumi_/1000.,"bRegTop.|top."));
      samples.push_back(MySample("Data, noBReg", "Run2012", kSigVar, kGreen-5, 9, lumi_/1000.*9./4.,"bRegTop.|top."));
	}

    
    //Scale before/after b-regression
    if(plotSelectedForPlotting.find("ScaleTests")!=plotSelectedForPlotting.end()){
      samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
      samples.push_back(MySample("t#bar{t}, Z2*, noBReg", "Summer12_TTJets1725_1.00", kSigVar, kCyan+1, 1, lumi_/1000.,"bRegTop.|top."));
      samples.push_back(MySample("t#bar{t}, Sc-Up", "Summer12_TTJets1725_scaleup", kSigVar, kMagenta+1, 9, lumi_/1000.*9./4.));
      samples.push_back(MySample("t#bar{t}, Sc-Up, noBReg", "Summer12_TTJets1725_scaleup", kSigVar, kBlue+1, 9, lumi_/1000.*9./4.,"bRegTop.|top."));
      samples.push_back(MySample("t#bar{t}, Sc-Dn", "Summer12_TTJets1725_scaledown", kSigVar, kGreen+1, 8, lumi_/1000.*9./4.));
      samples.push_back(MySample("t#bar{t}, Sc-Dn, noBReg", "Summer12_TTJets1725_scaledown", kSigVar, kAzure-2, 8, lumi_/1000.*9./4.,"bRegTop.|top."));
	}
	  
    //Matching before/after b-regression
	if(plotSelectedForPlotting.find("MatchingTests")!=plotSelectedForPlotting.end()){
	  samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
	  samples.push_back(MySample("t#bar{t}, Z2*, noBReg", "Summer12_TTJets1725_1.00", kSigVar, kCyan+1, 1, lumi_/1000.,"bRegTop.|top."));
	  samples.push_back(MySample("t#bar{t}, M-Up", "Summer12_TTJets1725_matchingup", kSigVar, kMagenta+1, 9, lumi_/1000.*9./4.));
	  samples.push_back(MySample("t#bar{t}, M-Up, noBReg", "Summer12_TTJets1725_matchingup", kSigVar, kBlue+1, 9, lumi_/1000.*9./4.,"bRegTop.|top."));
	  samples.push_back(MySample("t#bar{t}, M-Dn", "Summer12_TTJets1725_matchingdown", kSigVar, kGreen+1, 8, lumi_/1000.*9./4.));
	  samples.push_back(MySample("t#bar{t}, M-Dn, noBReg", "Summer12_TTJets1725_matchingdown", kSigVar, kAzure-2, 8, lumi_/1000.*9./4.,"bRegTop.|top."));
	}

//    /* Signal variations (plaintests)*/
//    samples.push_back(MySample("t#bar{t}, Z2*", "Summer12_TTJets1725_1.00", kSigVar, kRed+1, 1, lumi_/1000.));
//    samples.push_back(MySample("t#bar{t}, P11", "Summer12_TTJets1725_MGDecays_P11", kSigVar, kMagenta+1, 9, lumi_/1000.*9./4.));
////    samples.push_back(MySample("t#bar{t}, Powheg+Pythia", "Summer12_Z2_S12_POWHEG_sig_1.00", kSigVar, kGreen+1, 9, lumi_/1000.));
////    samples.push_back(MySample("t#bar{t}, Powheg+Herwig", "Summer12_Z2_S12_POWHER_sig_1.00", kSigVar, kBlue+1, 7, lumi_/1000.));
////    samples.push_back(MySample("t#bar{t}, MC@NLO", "Summer12_TTJets1725_mcatnlo_herwig", kSigVar, kCyan+1, 1, lumi_/1000.));
    //*/

  }
  
  //
  // FILL HISTOGRAMS
  //

  {
    int bkgCounter = -1;
    int sigVarCounter = -1;
    int bkgVarCounter = -1;
    // Loop over all samples
    std::cout << "Adding files from folder: " << path_ << std::endl;
    for(MySample& sample : samples){
      if(sample.type == kBkg    ) ++bkgCounter;
      if(sample.type == kSigVar ) ++sigVarCounter;
      if(sample.type == kBkgVar ) ++bkgVarCounter;
      // Get sample
      TChain* chain; int nFiles = 0;
      if (channelID == Helper::kAllJets) {
        chain = new TChain("analyzeKinFit/eventTree");
        //if(sample.type == kBkg)
        //  nFiles = chain->Add((sample.file+std::string(".root")).c_str());
        //else
        nFiles = chain->Add((path_+sample.file+std::string(".root")).c_str());
      }
      else if(plotSelectedForPlotting.find("GBRTesting")!=plotSelectedForPlotting.end()){
    	  chain = new TChain("hmvacor");
    	  nFiles += chain->Add((path_+"fmvacorNew"+sample.file+std::string("BRegJet_genJetPt_BRegJet_jetPtCorr.root")).c_str());
    	  std::cout << "added hmvacor chain" << std::endl;
      }
      else if(plotSelectedForPlotting.find("CalibTreeMaker")!=plotSelectedForPlotting.end()){
    	  chain = new TChain("DiJetTree");
    	  nFiles += chain->Add((path_+sample.file+"/ak5PFCHS_*.root").c_str());
    	  std::cout << "added CalibTreeMaker kind of chain" << std::endl;
      }
      else if(plotSelectedForPlotting.find("eventTreeWithJEC")!=plotSelectedForPlotting.end()){
    	  chain = new TChain("eventTree");
    	  if (channelID != Helper::kElectronJets) nFiles += chain->Add((path_+sample.file+std::string("_muon/job_*.root")).c_str());
    	  if (channelID != Helper::kMuonJets) nFiles += chain->Add((path_+sample.file+std::string("_electron/job_*.root")).c_str());
      }
      else{
    	  chain = new TChain("analyzeHitFit/eventTree");
    	  if (channelID != Helper::kElectronJets) nFiles += chain->Add((path_+sample.file+std::string("_muon/job_*.root")).c_str());
    	 if (channelID != Helper::kMuonJets    ) nFiles += chain->Add((path_+sample.file+std::string("_electron/job_*.root")).c_str());
      }
      std::cout << "Adding " << nFiles << " files for " << sample.name << " (" << sample.file << "), type: " << sample.type << std::endl;

      // Soft-initialize histograms and add sample
      for(MyHistogram& hist : hists) {
        if(hist.Dimension() == -1) continue;
        hist.Init(chain, topBranchName_,sample.replaceVar);
        if     (sample.type == kData   ) hist.SetupData(&sample);//dataContainsMC bool will be set to true in case the sample name does not contain "Data"
        else if(sample.type == kSig    ) hist.AddSignal(&sample);
        else if(sample.type == kBkg    ) hist.AddBackground(&sample);
        else if(sample.type == kSigVar ) hist.AddSignalVariation(&sample);
        else if(sample.type == kBkgVar ) hist.AddBackgroundVariation(&sample);
      }

      //replaceVar used for replacements in selection

      // Initialize weight and selection formulas
      //TTreeFormula weight("weight",  po::GetOption<std::string>("weight").c_str(), chain);
      TTreeFormula weight("weight",  po::GetOptionReplaced("weight",sample.replaceVar).c_str(), chain);
      std::string tempSel(po::GetOptionReplaced("analysisConfig.selection",sample.replaceVar));
      TTreeFormula sel   ("sel"   ,  tempSel.c_str(), chain);
      TTreeFormula selCP ("selCP" , (po::GetOptionReplaced("analysisConfig.selection",sample.replaceVar)
				 +std::string(" & ")+po::GetOptionReplaced("analysisConfig.selectionCP",sample.replaceVar)).c_str(), chain);
      TTreeFormula selWP ("selWP" , (po::GetOptionReplaced("analysisConfig.selection",sample.replaceVar)
				 +std::string(" & ")+po::GetOptionReplaced("analysisConfig.selectionWP",sample.replaceVar)).c_str(), chain);
      TTreeFormula selUN ("selUN" , (po::GetOptionReplaced("analysisConfig.selection",sample.replaceVar)
                 +std::string(" & ")+po::GetOptionReplaced("analysisConfig.selectionUN",sample.replaceVar)).c_str(), chain);


      //  Loop over all events
      int    loopSize     = chain->GetEntries();
      double loopFraction = 1.;
      if (po::GetOption<bool>("test")) loopFraction = 0.01;
      for(int i = 0; i < loopSize*loopFraction; ++i){
    	  if(i%1000==0) std::cout<<"\r  Loop over permutations: "<<(int)((float)i*100./(float)(loopSize*loopFraction-1))<<"%";
        //if(i%1000==0)std::cout << "event " << (Long_t) i << std::endl;
        long entry = chain->LoadTree(i);
        if(entry  < 0) break;
        //if(i > 10000) break;
        if(entry == 0){
          for(auto& hist : hists) {
            for(auto& var : hist.varx) var->UpdateFormulaLeaves();
            for(auto& var : hist.vary) var->UpdateFormulaLeaves();
            if (hist.selection.size() > 0) hist.sel->UpdateFormulaLeaves();
            hist.histoweight->UpdateFormulaLeaves();
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
            for(auto& var : hist.varx) if(!var->GetNdata()) continue;
            if(hist.Dimension() == -1) continue;
            for(auto& var : hist.varx) if(var->GetNdata()<=j) continue;
            if(!hist.histoweight->GetNdata())continue;
            if(hist.Dimension() == 1){
              for(auto& var : hist.vary) if(!var->GetNdata()) continue;
            }
            else if(hist.Dimension() == 1){
              for(auto& var : hist.vary) if(var->GetNdata()<=j) continue;
            }
	    else if(hist.Dimension() == 2){
	      for(auto& var : hist.vary) if(!var->GetNdata()) continue;
	    }
            // Skip permutation if individual selection is not met
            if (hist.selection.size() > 0) {
              if(!hist.sel->GetNdata()) continue;
              if(!hist.sel->EvalInstance(j)) continue;
            }
	    // Skip permutations except the first one if "Events" is contained in the y-axis title
	    if (TString(hist.Data1D()->GetYaxis()->GetTitle()).Contains("Events") && j > 0) continue;
            // Fill according to sample and histogram type
            // Fill data
            if     (sample.type == kData && hist.Dimension() == 1) {
              for(auto& var : hist.varx) hist.Data1D()->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            }
            else if(sample.type == kData && hist.Dimension() == 2) {
              for(auto& var : hist.varx) for(auto& var2 : hist.vary) hist.Data2D()->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            }
            // Fill signal
            else if(sample.type == kSig && hist.Dimension() == 1) {
              for(auto& var : hist.varx) {
                if     (selCP.EvalInstance(j)) hist.Sig1D().at(0)->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                else if(selWP.EvalInstance(j)) hist.Sig1D().at(1)->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                else if(selUN.EvalInstance(j)) hist.Sig1D().at(2)->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
              }
            }
            else if(sample.type == kSig && hist.Dimension() == 2) {
              for(auto& var : hist.varx) {
                for(auto& var2 : hist.vary) {
                    if     (selCP.EvalInstance(j)) hist.Sig2D().at(0)->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                    else if(selWP.EvalInstance(j)) hist.Sig2D().at(1)->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                    else if(selUN.EvalInstance(j)) hist.Sig2D().at(2)->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                }
              }
            }
            // Fill background
            else if(sample.type == kBkg  && hist.Dimension() == 1) {
              for(auto& var : hist.varx) {
                hist.Bkg1D()[bkgCounter]->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
                // Add background also to histogram for signal variation
                // TODO: Background normalization in signal variations not correct for AllJets
                //for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
              }
            }
            else if(sample.type == kBkg  && hist.Dimension() == 2) {
              for(auto& var : hist.varx) for(auto& var2 : hist.vary) hist.Bkg2D()[bkgCounter]->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
              // Add background also to histogram for signal variation
              // TODO: Background normalization in signal variations not correct for AllJets
              for(TH2F* sigvar : hist.Sigvar2D()) for(auto& var : hist.varx) for(auto& var2 : hist.vary) sigvar->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            }
            // Fill signal variation
            else if(sample.type == kSigVar && hist.Dimension() == 1) for(auto& var : hist.varx) hist.Sigvar1D()[sigVarCounter]->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            else if(sample.type == kSigVar && hist.Dimension() == 2) for(auto& var : hist.varx) for(auto& var2 : hist.vary) hist.Sigvar2D()[sigVarCounter]->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            // Fill background variation
            else if(sample.type == kBkgVar && hist.Dimension() == 1) for(auto& var : hist.varx) hist.Bkgvar1D()[bkgVarCounter]->Fill(var->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
            else if(sample.type == kBkgVar && hist.Dimension() == 2) for(auto& var : hist.varx) for(auto& var2 : hist.vary) hist.Bkgvar2D()[bkgVarCounter]->Fill(var->EvalInstance(j), var2->EvalInstance(j), weight.EvalInstance(j)*hist.histoweight->EvalInstance(j)*sample.scale);
          }
        }
      }
      for(MyHistogram& hist : hists) {
        for(auto& var : hist.varx) var->Clear();
        for(auto& var : hist.vary) var->Clear();
      }
    }
  }
  
  std::cout << "ended filling histograms" << std::endl;
  
  // 
  // ADD BACKGROUND TO VARIATIONS
  //
  
  for(MyHistogram& hist : hists) {
    if(hist.Dimension() != 1) continue;
    
    for(TH1F* sigvar : hist.Sigvar1D()) {
      for(TH1F* bkg : hist.Bkg1D()) {
        sigvar->Add(bkg);
      }
    }
    
    // Add normalization uncertainties to histos
    // TODO: double-counting of norm uncertainties, normalizeToData impossible
    //int bkgCounter = -1;
    //for(MySample& sample : samples){
    //  if(sample.type == kBkg) {
    //    ++bkgCounter;
    //    TH1F* bkg = hist.Bkg1D()[bkgCounter];
    //    for(int i = 0; i < bkg->GetNbinsX()+2; ++i) bkg->SetBinError(i, sqrt(pow(bkg->GetBinError(i),2) + pow(bkg->GetBinContent(i)*sample.scaleunc,2)));
    //  }
    //  
    //  if(sample.type == kSig) {
    //    for(TH1F* sig : hist.Sig1D()) {
    //      for(int i = 0; i < sig->GetNbinsX()+2; ++i) sig->SetBinError(i, sqrt(pow(sig->GetBinError(i),2) + pow(sig->GetBinContent(i)*sample.scaleunc,2)));
    //    }
    //  }
    //}

    // Loop over all samples to add background norm variations
    int bkgCounter = -1;
    for(MySample& sample : samples){
      if(sample.type == kBkg) ++bkgCounter;
      
      if(sample.type == kSig) {
        hist.AddNormVariation(&sample, "nominal");
        for(TH1F* sig : hist.Sig1D()) hist.Sigvar1D().back()->Add(sig);
        for(TH1F* bkg : hist.Bkg1D()) hist.Sigvar1D().back()->Add(bkg);
      }
      
      if(sample.scaleunc == 0) continue;
      
      hist.AddNormVariation(&sample, "up");
      for(TH1F* sig : hist.Sig1D()) hist.Sigvar1D().back()->Add(sig);
      for(TH1F* bkg : hist.Bkg1D()) hist.Sigvar1D().back()->Add(bkg);
      
      if(sample.type == kSig) {
        for(TH1F* sig : hist.Sig1D()) hist.Sigvar1D().back()->Add(sig, sample.scaleunc);
      }
      if(sample.type == kBkg) {
        hist.Sigvar1D().back()->Add(hist.Bkg1D()[bkgCounter], sample.scaleunc);
      }
      
      hist.AddNormVariation(&sample, "down");
      for(TH1F* sig : hist.Sig1D()) hist.Sigvar1D().back()->Add(sig);
      for(TH1F* bkg : hist.Bkg1D()) hist.Sigvar1D().back()->Add(bkg);
      
      if(sample.type == kSig) {
        for(TH1F* sig : hist.Sig1D()) hist.Sigvar1D().back()->Add(sig, -sample.scaleunc);
      }
      if(sample.type == kBkg) {
        hist.Sigvar1D().back()->Add(hist.Bkg1D()[bkgCounter], -sample.scaleunc);
      }
    }
  }
  
  // YIELDS and NORMALIZATION
  
  bool firstHist = true;
  for(MyHistogram& hist : hists) {
    if(hist.Dimension() != 1) continue;
    
    if (hist.name == "combinationType") {
      std::cout << "PERMUTATION YIELD" << std::endl;
      firstHist = true;
    }

    // Show event yields for first histogram
    
    int bins = hist.Data1D()->GetNbinsX()+2;
    double integral;
    double error = 0;
    double errorD = 0;
    double integralD = hist.Data1D()->Integral(0,bins);
    for(int i = 0; i < bins; ++i) errorD += hist.Data1D()->GetBinError(i);
    double errorS = 0;
    double integralS = 0;
    double errorB = 0;
    double integralB = 0;
    
    for(TH1F* sig : hist.Sig1D()) {
      integralS += sig->Integral(0,bins);
      //error = 0;
      //for(int i = 0; i < bins; ++i) error += sig->GetBinError(i);
      //errorS = sqrt(pow(errorS,2)+pow(error,2));
    }
    
    for(TH1F* bkg : hist.Bkg1D()) {
      integralB += bkg->Integral(0,bins);
      //error = 0;
      //for(int i = 0; i < bins; ++i) error += bkg->GetBinError(i);
      //errorB = sqrt(pow(errorB,2)+pow(error,2));
    }
    
    double integralSB = 0;
    if (channelID == Helper::kAllJets) integralSB = integralD;
    else                               integralSB = integralS+integralB;
    
    if (firstHist){
      std::cout << "======" << std::endl;
      std::cout << "Yields" << std::endl;
      std::cout << "======" << std::endl;
      printf("Data: %10.1lf +/- %-10.1lf \n", integralD, errorD);
      
      // TODO: use sample scaleUnc here
      int bkgCounter = -1;
      for(MySample& sample : samples){
        if(sample.type == kSig) {
          for(TH1F* sig : hist.Sig1D()) {
            integral = sig->Integral(0,bins);
            error = 0;
            for(int i = 0; i < bins; ++i) error += sqrt(pow(sig->GetBinError(i),2) + pow(sig->GetBinContent(i)*sample.scaleunc,2));
            errorS = sqrt(pow(errorS,2)+pow(error,2));
            printf("  %-20s %10.1lf +/- %-10.1lf %5.1lf %5.1lf \n", sig->GetTitle(), integral, error, integral/integralSB*100., integral/integralS*100.);
          }
          printf("Signal: %10.1lf +/- %-10.1lf %5.1lf \n", integralS, errorS, integralS/integralSB*100);
        }
        
        if(sample.type == kBkg) {
          ++bkgCounter;
          TH1F* bkg = hist.Bkg1D()[bkgCounter];
          integral = bkg->Integral(0,bins);
          error = 0;
          for(int i = 0; i < bins; ++i) error += sqrt(pow(bkg->GetBinError(i),2) + pow(bkg->GetBinContent(i)*sample.scaleunc,2));
          errorB = sqrt(pow(errorB,2)+pow(error,2));
          printf("  %-20s %10.1lf +/- %-10.1lf %5.1lf \n", bkg->GetTitle(), integral, error, integral/integralSB*100);
        }
      }
      printf("Background: %10.1lf +/- %-10.1lf %5.1lf \n", integralB, errorB, integralB/integralSB*100);
      printf("MC Total: %10.1lf +/- %-10.1lf \n", integralS+integralB, sqrt(pow(errorS,2)+pow(errorB,2)));
      
      firstHist = false;
    }
    
    // Alljets: Normalize to data
    if (channelID == Helper::kAllJets) {
      //double fSig = po::GetOption<double>("templates.fSig");
      //for(TH1F* sig    : hist.Sig1D())    sig->Scale(    fSig *integralD/integralS);
      //for(TH1F* bkg    : hist.Bkg1D())    bkg->Scale((1.-fSig)*integralD/integralB);
      //for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Scale(integralD/sigvar->Integral(0,bins));
      
      for(TH1F* bkg    : hist.Bkg1D()) {
        // Remove background from uncertainty histos
        for (TH1F* sigvar : hist.Sigvar1D()) {
          sigvar->Add(bkg, -1.);
        }
        // Scale
        bkg->Scale((integralD-integralS)/integralB);
        //Add back
        for (TH1F* sigvar : hist.Sigvar1D()) {
          sigvar->Add(bkg);
        }
      }
      for(TH1F* bkgvar : hist.Bkgvar1D()) {
        // Scale
        double integralBvar = bkgvar->Integral(0,bins);
        bkgvar->Scale((integralD-integralS)/integralBvar);
        //Add signal
        for (TH1F* sig : hist.Sig1D()) {
          bkgvar->Add(sig);
        }
      }
    }
    // Lepton+jets: Normalize to data
    else {
      if (po::GetOption<bool>("analysisConfig.normalizeToData")) {
        std::cout << "Normalizing to data" << std::endl;
        for(TH1F* sig    : hist.Sig1D())    sig->Scale(integralD/(integralS+integralB));
        for(TH1F* bkg    : hist.Bkg1D())    bkg->Scale(integralD/(integralS+integralB));
        for(TH1F* sigvar : hist.Sigvar1D()) sigvar->Scale(integralD/sigvar->Integral(0,bins));
      }
      //*/
      // Double bin error due to correlations
      // (Many observables do not change for different neutrino solutions or b assignments)
      if (po::GetOption<bool>("analysisConfig.plotPermutations")) {
        for (int i = 0; i < hist.Data1D()->GetNbinsX(); ++i) {
          hist.Data1D()->SetBinError(i, hist.Data1D()->GetBinError(i) * 2);
        }
        for(TH1F* sig : hist.Sig1D()) {
          for (int i = 0; i < sig->GetNbinsX(); ++i) {
            sig->SetBinError(i, sig->GetBinError(i) * 2);
          }
        }
        for(TH1F* bkg : hist.Bkg1D()) {
          for (int i = 0; i < bkg->GetNbinsX(); ++i) {
            bkg->SetBinError(i, bkg->GetBinError(i) * 2);
          }
        }
      }
    }
  }
  
  //
  // CACLCULATE UNCERTAINTIES
  //
  
  for(MyHistogram& hist : hists) {
    if(hist.Dimension() != 1) continue;
    
    // Add all stacked histos to uncertainty histogram, reset their uncertainty
    for(TH1F* sig : hist.Sig1D()) {
      hist.Unc1D()->Add(sig);
      for(int i = 0; i < sig->GetNbinsX()+2; ++i) sig->SetBinError(i, 0);
    }
    for(TH1F* bkg : hist.Bkg1D()) {
      hist.Unc1D()->Add(bkg);
      for(int i = 0; i < bkg->GetNbinsX()+2; ++i) bkg->SetBinError(i, 0);
    }
    
    // Init uncertainty vectors
    std::vector<double> up(hist.Unc1D()->GetNbinsX()+2,0), down(hist.Unc1D()->GetNbinsX()+2,0);

    // Fill uncertainty vectors
    for(TH1F* sigvar : hist.Sigvar1D()) {
      for(int i = 0; i < hist.Unc1D()->GetNbinsX()+2; ++i) {
        if (sigvar->GetBinContent(i) > hist.Unc1D()->GetBinContent(i)) {
          up[i] = sqrt(pow(up[i], 2) + pow(sigvar->GetBinContent(i) - hist.Unc1D()->GetBinContent(i), 2));
        }
        else {
          down[i] = sqrt(pow(down[i], 2) + pow(sigvar->GetBinContent(i) - hist.Unc1D()->GetBinContent(i), 2));
        }
      }
    }
    
    for(TH1F* bkgvar : hist.Bkgvar1D()) {
      for(int i = 0; i < hist.Unc1D()->GetNbinsX()+2; ++i) {
        if (bkgvar->GetBinContent(i) > hist.Unc1D()->GetBinContent(i)) {
          up[i] = sqrt(pow(up[i], 2) + pow(bkgvar->GetBinContent(i) - hist.Unc1D()->GetBinContent(i), 2));
        }
        else {
          down[i] = sqrt(pow(down[i], 2) + pow(bkgvar->GetBinContent(i) - hist.Unc1D()->GetBinContent(i), 2));
        }
      }
    }
    
    // Set bin center and error for uncertainty band
    for(int i = 0; i < hist.Unc1D()->GetNbinsX()+2; ++i) {
      //hist.Unc1D()->SetBinContent(i, hist.Unc1D()->GetBinContent(i) + (up[i]-down[i])/2.);
      //hist.Unc1D()->SetBinError(i, sqrt(pow(hist.Unc1D()->GetBinError(i), 2) + pow((up[i]+down[i])/2., 2)));
      hist.Unc1D()->SetBinError(i, std::max(std::abs(up[i]),std::abs(down[i])));
    }
  }

  //
  // DRAW CONTROL PLOTS
  //

  gSystem->mkdir((std::string("plot/controlplots/").c_str()));

  // initialize file stream for logging
  std::ofstream logResultsFile;
  logResultsFile.open ((std::string("plot/controlplots/")+channel_+".txt").c_str());
  
  // initialize root file for optional plot saving
  TFile* outFile = new TFile((std::string("plot/controlplots/")+channel_+".root").c_str(),"RECREATE");

  for(MyHistogram& hist : hists){
    if(hist.Dimension() == -1) continue;
    
    // move exponent for y-axis
    //TGaxis::SetExponentOffset(-0.05, 0.01, "y");
    
    // Create directory
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2D")).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2DProfile")).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1D")).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DRatio")).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVar")).c_str());
    gSystem->mkdir((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVarRatio")).c_str());

    // Draw 2D plots
    if(hist.Dimension() == 2){
      std::cout << "doing 2D plots" << std::endl;

      TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
      canv->cd();
      canv->SetRightMargin(0.15);
      canv->SetLogx(hist.LogX());
      canv->SetLogy(hist.LogY());


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
    	  //std::cout << "collectAll2D.at(h_i)->GetName()+(TString)\"_pfx\" " << collectAll2D.at(h_i)->GetName()+(TString)"_pfx" << std::endl;
    	  collectAll2D.at(h_i)->Draw("colz");
    	  collectAll2DProfiles.at(h_i)->Draw("same");

    	  TLegend* legInd = new TLegend(0.25, 0.73, 0.55, 0.905);
    	  legInd->SetTextSize(0.04);
    	  legInd->SetFillStyle(0);
    	  legInd->SetBorderSize(0);
    	  legInd->AddEntry( collectAll2D.at(h_i), collectAll2D.at(h_i)->GetTitle(), "LP" );
    	  legInd->Draw();

    	  helper->DrawCMS(-1, -1, canv);
    	  gPad->RedrawAxis();
    	  //std::cout << HelperFunctions::cleanedName(collectAll2D.at(h_i)->GetTitle()) << std::endl;
    	  canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2D/")+channel_+outPath_+std::string("_")+std::string(hist.Data2D()->GetName())+std::string("_Ind2D_")+HelperFunctions::cleanedName(collectAll2D.at(h_i)->GetTitle())+std::string(".eps")).c_str(),"eps");
          canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2D/")+channel_+outPath_+std::string("_")+std::string(hist.Data2D()->GetName())+std::string("_Ind2D_")+HelperFunctions::cleanedName(collectAll2D.at(h_i)->GetTitle())+std::string(".png")).c_str(),"png");
      }

      HelperFunctions::setCommonYRange(collectAll2DProfiles,0.35);
      
      canv->SetRightMargin(0.05);

      for (size_t h_i=0;h_i<collectAll2DProfiles.size();++h_i){
        if(hist.plotRangeYMin!=-99)collectAll2DProfiles.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYMin, hist.plotRangeYMax);
        if(h_i==0)collectAll2DProfiles.at(h_i)->Draw();
        else collectAll2DProfiles.at(h_i)->Draw("SAME");
      }

      TLegend* leg0 = new TLegend(0.25, 0.73, 0.55, 0.905);
      leg0->SetTextSize(0.04);
      leg0->SetFillStyle(0);
      leg0->SetBorderSize(0);
      for(TH2F* sig : hist.Sig2D()){
        if(!(channelID == Helper::kAllJets && TString(sig->GetTitle()).Contains("unmatched"))) leg0->AddEntry( sig, sig->GetTitle(), "LP" );
      }
      leg0->Draw();

      TLegend *leg1 = new TLegend(0.6, 0.73, 0.9, 0.905);
      leg1->SetTextSize(0.03);
      leg1->SetFillStyle(0);
      leg1->SetBorderSize(0);
      for(TH2F* bkg : hist.Bkg2D()) leg1->AddEntry( bkg, bkg->GetTitle(), "LP" );
      leg1->AddEntry( hist.Data2D(), hist.Data2D()->GetTitle(), "LP" );
      leg1->Draw();

      helper->DrawCMS(-1, -1, canv);
      gPad->RedrawAxis();


      canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2D/")+channel_+outPath_+std::string("_")+std::string(hist.Data2D()->GetName())+std::string(".eps")).c_str(),"eps");
      canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2D/")+channel_+outPath_+std::string("_")+std::string(hist.Data2D()->GetName())+std::string(".png")).c_str(),"png");


      collectAll2D.clear();
      for(TH1* prof : collectAll2DProfiles) delete prof;//delete profiles
      collectAll2DProfiles.clear();

      // Draw signal variation plots only if there are signal variations
      if(hist.Sigvar2D().size()){
        collectAll2D.push_back(hist.Data2D());
        for(TH2F* sigvar : hist.Sigvar2D()) collectAll2D.push_back(sigvar);

        for (size_t h_i=0;h_i<collectAll2D.size();++h_i){
          //std::cout << " title " << collectAll2D.at(h_i)->GetTitle() << " " << collectAll2D.at(h_i)->GetFillColor() << std::endl;
          collectAll2D.at(h_i)->SetMarkerStyle(20+h_i);
          //collectAll2D.at(h_i)->SetLineColor(collectAll2D.at(h_i)->GetFillColor());
          collectAll2D.at(h_i)->SetMarkerColor(collectAll2D.at(h_i)->GetLineColor());
          if(h_i==0){
            collectAll2D.at(h_i)->SetMarkerColor(1);
            collectAll2D.at(h_i)->SetLineColor(1);
          }
          collectAll2DProfiles.push_back(collectAll2D.at(h_i)->ProfileX(collectAll2D.at(h_i)->GetName()+(TString)"_pfx")->ProjectionX()); //convert TProfile to TH1 with ProjectionX()
          collectAll2DProfiles.at(h_i)->SetLineColor(collectAll2D.at(h_i)->GetLineColor());
          collectAll2DProfiles.at(h_i)->SetMarkerColor(collectAll2D.at(h_i)->GetLineColor());
          collectAll2DProfiles.at(h_i)->SetMarkerStyle(collectAll2D.at(h_i)->GetMarkerStyle());
	  
          collectAll2DProfiles.at(h_i)->GetYaxis()->SetTitle(collectAll2D.at(h_i)->GetYaxis()->GetTitle());
          if(hist.ExportSigVarToRoot()){
            //std::cout << "this should go into a root file" << std::endl;
            outFile->WriteTObject(collectAll2DProfiles.back());
          }

        }

        HelperFunctions::setCommonYRange(collectAll2DProfiles,0.35);

        for (size_t h_i=0;h_i<collectAll2DProfiles.size();++h_i){
          if(hist.plotRangeYMin!=-99)collectAll2DProfiles.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYMin, hist.plotRangeYMax);
          if(h_i==0){
        	hist.DataContainsMC()==false ? collectAll2DProfiles.at(h_i)->Draw("p") : collectAll2DProfiles.at(h_i)->Draw("hist");
          }
          else {
            collectAll2DProfiles.at(h_i)->Draw("hist SAME");
          }
          if(hist.DataContainsMC()==true)collectAll2DProfiles.at(h_i)->Draw("p SAME"); // plot with markers in case there is only MC
        }

        leg1->Clear();
        for(TH2F* sigvar : hist.Sigvar2D()) leg1->AddEntry( sigvar, sigvar->GetTitle(), "LP" );
        leg1->AddEntry( hist.Data2D(), hist.Data2D()->GetTitle(), hist.DataContainsMC()==false ? "P" : "LP");
        leg1->Draw();

        helper->DrawCMS(-1, -1, canv);

        gPad->RedrawAxis();

        canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2DProfile/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.eps")).c_str(),"eps");
        canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2DProfile/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar.png")).c_str(),"png");

        TCanvas* canvWRatio = new TCanvas("cControlPlotsWRatio", "cControlPlotsWRatio", 600, 600);
        canvWRatio->Range(0,0,1,1);
        canvWRatio->cd();
        canvWRatio->SetLogx(hist.LogX());
        canvWRatio->SetLogy(hist.LogY());

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,1,10); //bottom
        pad2->SetFillStyle(4100);
        pad2->SetTopMargin(0.71);
        pad2->Draw();
        pad2->cd();
        pad2->SetGridy();
        pad2->SetLogx(hist.LogX());

        canvWRatio->cd();

        TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1,10); //top
        pad1->SetFillStyle(4100);
        pad1->SetBottomMargin(0.3);
        //pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogx(hist.LogX());
        pad1->SetLogy(hist.LogY());

        double oldLabelSize = hist.Data2D()->GetLabelSize();
        double oldTitleSize = hist.Data2D()->GetTitleSize();
        hist.Data2D()->SetLabelSize(0);
        hist.Data2D()->SetTitleSize(0);

        TH1* data2DProfile = (TH1*) hist.Data2D()->ProfileX(hist.Data2D()->GetName()+(TString)"_pfx")->ProjectionX();
        data2DProfile->SetLineColor(hist.Data2D()->GetLineColor());
        data2DProfile->GetYaxis()->SetTitle(hist.Data2D()->GetYaxis()->GetTitle());

        std::vector<TH1*> collectRatios;
        for(TH1* sigvar : collectAll2DProfiles){
          collectRatios.push_back(HelperFunctions::createRatioPlot(sigvar, data2DProfile, (std::string)"MC/"+hist.Data2D()->GetTitle()));
          collectRatios.back()->UseCurrentStyle();
          //std::cout << collectRatios.back()->GetName() << std::endl;
          collectRatios.back()->SetMarkerColor(sigvar->GetLineColor());
          collectRatios.back()->SetLineColor(sigvar->GetLineColor());
          collectRatios.back()->SetLineStyle(sigvar->GetLineStyle());
        }
        pad2->cd();
        pad2->Draw();
        //std::cout << "collectRatios.size() " << collectRatios.size() << std::endl;
        for (size_t h_i=0;h_i<collectRatios.size();++h_i){
	  collectRatios.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYRatioMin,hist.plotRangeYRatioMax); //needs to be set properly for all histos for Draw->("hist same") to work properly... slightly unexpected behavior
	  if(h_i==0){
	    collectRatios.at(h_i)->Draw("hist");
	    collectRatios.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYRatioMin,hist.plotRangeYRatioMax);
	    collectRatios.at(h_i)->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
	    collectRatios.at(h_i)->GetYaxis()->SetNdivisions(205);
	    collectRatios.at(h_i)->GetYaxis()->CenterTitle();
	    collectRatios.at(h_i)->SetMarkerSize(0);
	    collectRatios.at(h_i)->SetFillColor(collectRatios.at(h_i)->GetLineColor());
	    collectRatios.at(h_i)->SetFillStyle(3254);
	    collectRatios.at(h_i)->DrawClone("e2 same");
	    collectRatios.at(h_i)->SetFillStyle(0);
	  }
	  else collectRatios.at(h_i)->Draw("HIST SAME ");
        }
        pad1->cd();
        data2DProfile->GetYaxis()->UnZoom();
        data2DProfile->GetYaxis()->SetRangeUser(data2DProfile->GetMinimum()*0.9, data2DProfile->GetMaximum()*1.15);
        if(hist.plotRangeYMin!=-99)data2DProfile->GetYaxis()->SetRangeUser(hist.plotRangeYMin, hist.plotRangeYMax);

        hist.DataContainsMC()==false ? data2DProfile->Draw("p") : data2DProfile->Draw("hist");
        for(TH1* sigvar : collectAll2DProfiles) sigvar->Draw("hist same");
        for(TH1* sigvar : collectAll2DProfiles) sigvar->Draw("p same");
        hist.DataContainsMC()==false ? data2DProfile->Draw("p same") : data2DProfile->Draw("hist same");
        leg1->Draw();

        canvWRatio->cd();
        helper->DrawCMS(-1, -1, canvWRatio);
        gPad->RedrawAxis();

        canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2DProfile/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar_Ratio.eps")).c_str(),"eps");
        canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/2DProfile/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_sigvar_Ratio.png")).c_str(),"png");

        hist.Data2D()->SetLabelSize(oldLabelSize);
        hist.Data2D()->SetTitleSize(oldTitleSize);

        delete data2DProfile;
        for(TH1* sigvar : collectAll2DProfiles) delete sigvar;        
      }//end if 2D signal variation plots
    }
    else if(hist.Dimension() == 1){

      // Draw control plot
      TCanvas* canv = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
      canv->cd();
      canv->SetLogx(hist.LogX());
      canv->SetLogy(hist.LogY());
      hist.Data1D()->GetYaxis()->SetRangeUser(1, hist.Data1D()->GetMaximum()*1.6);
      if(hist.LogY())hist.Data1D()->GetYaxis()->SetRangeUser(2, 2 *pow(hist.Data1D()->GetMaximum()/2,1/0.7));
      hist.DataContainsMC()==false ? hist.Data1D()->Draw("p") : hist.Data1D()->Draw("hist");
      THStack* stack = new THStack("stack", "");


      //unfortunately PYTHON-style for loop breaks when doing boost::adaptors::reverse(hist.Bkg1D()) directly (pointers get mixed up)
      std::vector <TH1F*> tempHistBkg1D = hist.Bkg1D();
      for(TH1F* bkg : boost::adaptors::reverse(tempHistBkg1D)) stack->Add(bkg);
      for(TH1* sig : hist.Sig1D()) stack->Add(sig);

      //for(int i = 0; i< hist.Data1D()->GetNbinsX(); i++){
      //  std::cout << ((TH1*)stack->GetStack()->Last())->GetBinContent(i) << " error" << ((TH1*)stack->GetStack()->Last())->GetBinError(i) << " data" << hist.Data1D()->GetBinContent(i)  << "; " << std::endl;
      //}
      //hs->GetStack()->Last()->Draw();
      if(hist.PlotStackNorm()){ //work in progress... normalize stack plot bin-by-bin
        TH1* sumstack = (TH1*) stack->GetStack()->Last()->Clone();
        stack = new THStack("stackNORM", "");
        TH1* dataTemp = (TH1*) hist.Data1D()->Clone();

        hist.setData1D((TH1F*) HelperFunctions::createRatioPlot(dataTemp, dataTemp, (std::string)"fraction "+hist.Data1D()->GetYaxis()->GetTitle()));
        hist.DataContainsMC()==false ? hist.Data1D()->Draw("p") : hist.Data1D()->Draw("hist");
        hist.Data1D()->GetYaxis()->SetRangeUser(0.01, 1.6);
        hist.DataContainsMC()==false ? hist.Data1D()->Draw("p") : hist.Data1D()->Draw("hist");

        std::vector <TH1F*> TEMPHistBkg1D = hist.Bkg1D();
        for(TH1F* bkg : boost::adaptors::reverse(TEMPHistBkg1D)) stack->Add(HelperFunctions::createRatioPlot(bkg,sumstack,bkg->GetTitle()+(std::string)"_norm"));
        std::vector <TH1F*> TEMPHistSig1D = hist.Sig1D();
        if (channelID == Helper::kAllJets) {
          std::vector <TH1F*> TEMPHistSig1D = hist.Sig1D();
          for(TH1* sig : boost::adaptors::reverse(TEMPHistSig1D)) stack->Add(HelperFunctions::createRatioPlot(sig,sumstack,sig->GetTitle()+(std::string)"_norm"));
        }
        else for(TH1* sig : hist.Sig1D()) stack->Add(HelperFunctions::createRatioPlot(sig,sumstack,sig->GetTitle()+(std::string)"_norm"));

        std::cout << "trying to normalize stacks" << std::endl;
      }	
      std::cout << stack->GetName() << std::endl;
      stack->Draw("hist same");
      hist.DataContainsMC()==false ? hist.Data1D()->Draw("p same") : hist.Data1D()->Draw("hist same");

      // LEGEND
      int nSamples = 1; // data
      // signal
      if (po::GetOption<bool>("analysisConfig.plotPermutations")) {
        if(channelID == Helper::kAllJets) nSamples += 2;
        else                              nSamples += 3;
      }
      else nSamples += 1;
      // background
      nSamples += hist.Bkg1D().size();
      
      TLegend* leg0 = new TLegend(0.25, 0.73, 0.55, 0.905);
      leg0->SetTextSize(0.04);
      leg0->SetFillStyle(0);
      leg0->SetBorderSize(0);
      
      TLegend *leg1 = new TLegend(0.6, 0.73, 0.9, 0.905);
      leg1->SetTextSize(0.04);
      leg1->SetFillStyle(0);
      leg1->SetBorderSize(0);      
      
      if (po::GetOption<bool>("analysisConfig.plotPermutations")) {
        for(TH1F* sig : hist.Sig1D()) {
          if(!(channelID == Helper::kAllJets && TString(sig->GetTitle()).Contains("unmatched"))) leg0->AddEntry( sig, sig->GetTitle(), "F" );
        }
	if (channelID == Helper::kAllJets) {
	  for(TH1F* bkg : hist.Bkg1D()) leg1->AddEntry( bkg, bkg->GetTitle(), "F" );
	  hist.DataContainsMC()==false ? leg1->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "P" ) : leg1->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "L" );
	} else {
	  hist.DataContainsMC()==false ? leg0->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "P" ) : leg0->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "L" );
	  for(TH1F* bkg : hist.Bkg1D()) leg1->AddEntry( bkg, bkg->GetTitle(), "F" );
        }
      }
      else {
        TLegend *currentLeg = leg0;
        int iSample = 1;
        currentLeg->AddEntry( hist.Sig1D()[0], "t#bar{t}", "F" );
        for(TH1F* bkg : hist.Bkg1D()) {
          ++iSample;
          if (iSample > nSamples/2.) currentLeg = leg1;
          currentLeg->AddEntry( bkg, bkg->GetTitle(), "F" );
        }
        hist.DataContainsMC()==false ? currentLeg->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "P" ) : currentLeg->AddEntry( hist.Data1D(), hist.Data1D()->GetTitle(), "L" );
      }
      
      leg0->Draw();
      leg1->Draw();

      gPad->RedrawAxis();
      helper->DrawCMS(-1, -1, canv);

      canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1D/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string(".eps")).c_str(),"eps");
      canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1D/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string(".png")).c_str(),"png");


      //draw together with ratio underneath
      TH1* ratioToTHStack = HelperFunctions::createRatioPlot((TH1*)hist.Data1D(), ((TH1*)(stack->GetStack()->Last())), hist.Data1D()->GetTitle()+(std::string)"/MC");
      TH1* ratioUncTHStack = HelperFunctions::createRatioPlot((TH1*)hist.Unc1D(), ((TH1*)(stack->GetStack()->Last())), hist.Unc1D()->GetTitle()+(std::string)"/MC");
      TCanvas* canvWRatio = new TCanvas("cControlPlotsWRatio", "cControlPlotsWRatio", 600, 600);
      canvWRatio->Range(0,0,1,1);
      canvWRatio->cd();
      canvWRatio->SetLogx(hist.LogX());
      canvWRatio->SetLogy(hist.LogY());

      TPad *pad2 = new TPad("pad2","pad2",0,0,1,1,10); //bottom
      pad2->SetFillStyle(4100);
      pad2->SetFillStyle(4000);
      pad2->SetTopMargin(0.71);
      pad2->Draw();
      pad2->cd();
      pad2->SetGridy();
      pad2->SetLogx(hist.LogX());
      ratioToTHStack->Draw();
      ratioToTHStack->GetXaxis()->SetNdivisions(506);

      //std::cout << "hist.plotRangeYRatioMin" << hist.plotRangeYRatioMin << " hist.plotRangeYRatioMax " << hist.plotRangeYRatioMax << std::endl;
      ratioToTHStack->GetYaxis()->SetRangeUser(hist.plotRangeYRatioMin,hist.plotRangeYRatioMax);
      ratioToTHStack->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
      ratioToTHStack->GetYaxis()->SetNdivisions(205);
      ratioToTHStack->GetYaxis()->CenterTitle();
      
      ratioUncTHStack->Draw("e2 same");

      canvWRatio->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1,10); //top
      //pad1->SetFillStyle(4100);
      pad1->SetFillStyle(4000);
      pad1->SetBottomMargin(0.3);
      //pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      pad1->Draw();
      pad1->cd();
      pad1->SetLogx(hist.LogX());
      pad1->SetLogy(hist.LogY());
      
      double oldLabelSize = hist.Data1D()->GetLabelSize();
      double oldTitleSize = hist.Data1D()->GetTitleSize();
      hist.Data1D()->SetLabelSize(0);
      hist.Data1D()->SetTitleSize(0);

      if(hist.LogY())hist.Data1D()->GetYaxis()->SetRangeUser(2, 2 *pow(hist.Data1D()->GetMaximum()/2,1/0.7));

      
      hist.DataContainsMC()==false ? hist.Data1D()->Draw("p") : hist.Data1D()->Draw("hist");
      //hist.Data1D()->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
      hist.Data1D()->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y"));
      stack->Draw("hist same");
      hist.Unc1D()->Draw("e2 same");
      hist.DataContainsMC()==false ? hist.Data1D()->Draw("p same") : hist.Data1D()->Draw("hist same");
      //leg0->SetY1NDC(0.675);
      leg0->Draw();
      //leg1->SetY1NDC(0.675);
      leg1->Draw();
      
      gPad->RedrawAxis();

      canvWRatio->cd();
      
      helper->DrawCMS(-1, -1, canvWRatio);
      helper->DrawLabel(hist.extraLabel, 0.6, 0.63, 0.9);

      canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DRatio/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_Ratio.eps")).c_str(),"eps");
      canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DRatio/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_Ratio.png")).c_str(),"png");
      canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DRatio/")+channel_+outPath_+std::string("_")+std::string(hist.Data1D()->GetName())+std::string("_Ratio.root")).c_str(),"root");

      hist.Data1D()->SetLabelSize(oldLabelSize);
      hist.Data1D()->SetTitleSize(oldTitleSize);

      // Draw signal variation plot only if there are signal variations
      if(hist.Sigvar1D().size()){
        std::cout << "plotting sigvar" << std::endl;
        TH1F* hist1DData = hist.Data1D();
        std::vector<TH1F*> hist1DSigVars;
        std::vector<double> fSigVars;

        double sigIntegral = 0;
        for(TH1F* sig : hist.Sig1D()){
          sigIntegral += sig->Integral(0,sig->GetNbinsX()+2);
        }
        double fSig = sigIntegral/hist1DData->Integral(0,hist1DData->GetNbinsX()+2);

        double maxForPlot = hist1DData->GetMaximum();
          for(TH1F* sigvar : hist.Sigvar1D()){
            if(sigvar->GetMaximum() > maxForPlot) maxForPlot = sigvar->GetMaximum();
            hist1DSigVars.push_back((TH1F*)sigvar->Clone(sigvar->GetName()));
        }

        //if(channelID == Helper::kAllJets){
        //  hist1DData = (TH1F*)hist.Data1D()->Clone(hist.Data1D()->GetName());
        //  hist1DData->Add(hist.Bkg1D().at(0),-1.);
        //}
        //else{
        //  hist1DData = hist.Data1D();
        //}

        if(channelID == Helper::kAllJets){
          for(TH1F* sigvar : hist1DSigVars){
            double varIntegral = sigvar->Integral(0,sigvar->GetNbinsX()+2);
            double fSigVar = varIntegral/sigIntegral;
            fSigVars.push_back(fSigVar);
            sigvar->Add(hist.Bkg1D().at(0),(hist1DData->Integral(0,hist1DData->GetNbinsX()+2)-varIntegral)/hist.Bkg1D().at(0)->Integral(0,hist.Bkg1D().at(0)->GetNbinsX()+2));
            //std::cout << "data: "<< hist1DData->Integral(0,hist1DData->GetNbinsX()+2) << ", sigvar: " << varIntegral << ", bkg: " << hist.Bkg1D().at(0)->Integral(0,hist.Bkg1D().at(0)->GetNbinsX()+2) << ", fsig: " << fSig << ", fsigvar: " << fSigVar << ", sigInt: " << sigIntegral << ", after: " << sigvar->Integral(0,sigvar->GetNbinsX()+2) << std::endl;
          }
        }

        canv->cd();
        canv->Clear();

        hist1DData->GetYaxis()->UnZoom();
        if(hist.DataContainsMC()==false) hist1DData->Draw("p");
        else {
          hist1DData->Draw("hist");
          hist1DData->SetFillStyle(3445);
          hist1DData->SetFillColor(hist1DData->GetLineColor());
        }

        hist1DData->GetYaxis()->SetRangeUser(1, 1.5*maxForPlot);
        if(hist.LogY())hist.Data1D()->GetYaxis()->SetRangeUser(2, 2 *pow(hist.Data1D()->GetMaximum()/2,1/0.7));


        for(TH1F* sigvar : hist1DSigVars) sigvar->DrawClone("hist same");
        for(TH1F* sigvar : hist1DSigVars){
          sigvar->SetMarkerSize(0);
          sigvar->SetFillColor(sigvar->GetLineColor());
          sigvar->SetFillStyle(3254);
          sigvar->DrawClone("e2 same");
          sigvar->SetFillStyle(0);
        }
        hist.DataContainsMC()==false ? hist1DData->Draw("p same") : hist1DData->Draw("hist same");

        leg1->Clear();
        if(hist.legendHeader!="") leg1->SetHeader(hist.legendHeader.c_str()); 
        leg1->SetY1NDC(0.7);
        leg1->SetX1NDC(0.25);
        leg1->SetX2NDC(0.9);
        leg1->SetMargin(0.125);
        //leg1->SetX1NDC(0.63);
        //leg1->SetX2NDC(0.9);
        //leg1->SetMargin(0.3);
        gPad->Update();
        gPad->Modified();

        if(hist.FitGaussToCore()){
          TF1* gauss;
          double width, widthErr, rms, rmsErr;
          char buffer [50];
          std::string logResultsCSV = hist1DData->GetName()+std::string(", ") + HelperFunctions::cleanedName(hist1DData->GetXaxis()->GetTitle())+std::string(", ");
          for(TH1F* sigvar : hist1DSigVars){
            HelperFunctions::fitCoreWidth(sigvar, 1.5, gauss, width, widthErr, rms, rmsErr);
            sprintf(buffer,"%s, #sigma/#mu=%5.3f,#mu=%5.3f,m=%5.3f",sigvar->GetTitle(),width/gauss->GetParameter(1),gauss->GetParameter(1),sigvar->GetMean());
            leg1->AddEntry( sigvar, buffer, "L" );
            sprintf(buffer,"%s, %5.3f, %5.3f, %5.3f, ",HelperFunctions::cleanedName(sigvar->GetTitle()).c_str(),width/gauss->GetParameter(1),gauss->GetParameter(1),sigvar->GetMean());
            logResultsCSV+=buffer;
          }
          HelperFunctions::fitCoreWidth(hist1DData, 1.5, gauss, width, widthErr, rms, rmsErr);
          sprintf(buffer,"%s, #sigma/#mu=%5.3f,#mu=%5.3f,m=%5.3f",hist1DData->GetTitle(),width/gauss->GetParameter(1),gauss->GetParameter(1),hist1DData->GetMean());
          hist.DataContainsMC()==false ? leg1->AddEntry( hist1DData, buffer, "P" ) : leg1->AddEntry( hist1DData, buffer, "L" );
          sprintf(buffer,"%s, %5.3f, %5.3f, %5.3f \n",HelperFunctions::cleanedName(hist1DData->GetTitle()).c_str(),width/gauss->GetParameter(1),gauss->GetParameter(1),hist1DData->GetMean());
          logResultsCSV+=buffer;
          logResultsFile << logResultsCSV.c_str();
          delete gauss;

        }
        else {
          if(channelID == Helper::kAllJets) {
            for(unsigned int i = 0; i < hist1DSigVars.size(); ++i) {
              char temp[99];
              sprintf (temp, "%s + Bkg, f_{sig} = %4.1f%%", hist1DSigVars[i]->GetTitle(), 100.*fSig*fSigVars.at(i));
              //sprintf (temp, "%s", hist1DSigVars[i]->GetTitle());
              leg1->AddEntry( hist1DSigVars[i], temp, "L" );
            }
          }
          else for(TH1F* sigvar : hist1DSigVars) leg1->AddEntry( sigvar, sigvar->GetTitle(), "L" );
                hist.DataContainsMC()==false ? leg1->AddEntry( hist1DData, hist1DData->GetTitle(), "P" ) : leg1->AddEntry( hist1DData, hist1DData->GetTitle(), "L" );
        }
        outFile->WriteTObject(hist1DData);
        for(TH1F* sigvar : hist1DSigVars) outFile->WriteTObject(sigvar);
        leg1->Draw();

        helper->DrawCMS(-1, -1, canv);

        hist1DData->Draw("axissame");

        canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVar/")+channel_+outPath_+std::string("_")+std::string(hist1DData->GetName())+std::string("_sigvar.eps")).c_str(),"eps");
        canv->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVar/")+channel_+outPath_+std::string("_")+std::string(hist1DData->GetName())+std::string("_sigvar.png")).c_str(),"png");

        //ratio plots
        canvWRatio->cd();
        //canvWRatio->Clear();
        std::vector<TH1*> collectRatios;
        collectRatios.push_back(HelperFunctions::createRatioPlot(hist1DData, hist1DData, (std::string)"MC/"+hist.Data2D()->GetTitle()));
        for(TH1F* sigvar : hist1DSigVars){
          collectRatios.push_back(HelperFunctions::createRatioPlot(((TH1 *)(sigvar)), (TH1 *)hist1DData, (std::string)"MC/"+hist1DData->GetTitle()));
          collectRatios.back()->SetMarkerColor(sigvar->GetLineColor());
          collectRatios.back()->SetLineColor  (sigvar->GetLineColor());
          collectRatios.back()->SetLineStyle  (sigvar->GetLineStyle());
          collectRatios.back()->SetLineWidth  (sigvar->GetLineWidth());
        }

        pad2->cd();
        pad2->Draw();
        for (size_t h_i=0;h_i<collectRatios.size();++h_i){
          collectRatios.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYRatioMin,hist.plotRangeYRatioMax);
          if(h_i==0){
            collectRatios.at(h_i)->Draw("hist");
            collectRatios.at(h_i)->GetYaxis()->SetRangeUser(hist.plotRangeYRatioMin,hist.plotRangeYRatioMax);
            collectRatios.at(h_i)->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
            collectRatios.at(h_i)->GetYaxis()->SetNdivisions(205);
            collectRatios.at(h_i)->GetYaxis()->CenterTitle();
            collectRatios.at(h_i)->SetMarkerSize(0);
            collectRatios.at(h_i)->SetFillColor(collectRatios.at(h_i)->GetLineColor());
            collectRatios.at(h_i)->SetFillStyle(3254);
            collectRatios.at(h_i)->DrawClone("e2 same");
            collectRatios.at(h_i)->SetFillStyle(0);
          }
          else collectRatios.at(h_i)->DrawClone("hist same");
        }
        collectRatios.at(0)->Draw("axissame");
        pad1->Draw();
        pad1->cd();
        double oldLabelSize = hist1DData->GetLabelSize();
        double oldTitleSize = hist1DData->GetTitleSize();
        hist1DData->SetLabelSize(0);
        hist1DData->SetTitleSize(0);
        hist1DData->GetYaxis()->UnZoom();
        hist1DData->GetYaxis()->SetRangeUser(1, maxForPlot*1.75); //was 0.75 for thesis plot of nVertex
        if(hist.LogY())hist.Data1D()->GetYaxis()->SetRangeUser(2, 2 *pow(hist.Data1D()->GetMaximum()/2,1/0.7));
        hist.DataContainsMC()==false ? hist1DData->Draw("p") : hist1DData->Draw("hist");
        for(TH1F* sigvar : hist1DSigVars) sigvar->Draw("hist same");
        hist.DataContainsMC()==false ? hist1DData->Draw("p same") : hist1DData->Draw("hist same");
	hist1DData->Draw("axissame");
        leg1->Draw();

        canvWRatio->cd();
        helper->DrawCMS(-1, -1, canv);
        gPad->RedrawAxis();

        canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVarRatio/")+channel_+outPath_+std::string("_")+std::string(hist1DData->GetName())+std::string("_sigvar_Ratio.eps")).c_str(),"eps");
        canvWRatio->Print((std::string("plot/controlplots/")+channel_+outPath_+std::string("/1DVarRatio/")+channel_+outPath_+std::string("_")+std::string(hist1DData->GetName())+std::string("_sigvar_Ratio.png")).c_str(),"png");
        hist1DData->SetLabelSize(oldLabelSize);
	hist1DData->SetTitleSize(oldTitleSize);
      }//end if 1D signal variation plots
    }//end 1D plots
    else{
      std::cout << "The Histogram *" << hist.Data1D()->GetName() << "* has an invalid dimension *" << hist.Dimension() << "*! Supported are only 1 and 2!" << std::endl;
    }
  }
  logResultsFile.close();
  outFile->Close();
  
  std::_Exit(0);
}

void TopMassControlPlots::MyHistogram::ConfigureExtraOptions(bool SetFitGaussToCore, std::string SetCustomHistoweight, bool ExportSigVarToRoot, std::string LegendHeader, bool useLogXBins, bool setLogYOnPads, bool useLogYBins)
{
	fitGaussToCore = SetFitGaussToCore;
	CustomHistoWeight = SetCustomHistoweight;
	exportSigVarToRoot = ExportSigVarToRoot;
	legendHeader = LegendHeader;

	unsigned int nBinsX = data->GetNbinsX();
	double xmin = data->GetXaxis()->GetXmin();
	double xmax = data->GetXaxis()->GetXmax();
	std::vector<double> binEdgesX(nBinsX+1);
	if(useLogXBins)logX = true;
	if(setLogYOnPads)logY = true;
	HelperFunctions::equidistLogBins(binEdgesX, xmin, xmax, useLogXBins);

	if(Dimension()==1){
    	if(useLogXBins)data->SetBins(nBinsX,&binEdgesX[0]);
    	if(useLogYBins){
    			std::cout << "Error: It does not make sense to request log y bins for 1D histograms" << std::endl;
//    			std::flush;
    			std::_Exit(99);
    	}
	}
	else if(Dimension()==2){
//		for(double binlowedge : binEdgesX) std::cout <<", " << binlowedge;
//		std::cout << std::endl;
		unsigned int nBinsY = data->GetNbinsY();
		double ymin = data->GetYaxis()->GetXmin();
		double ymax = data->GetYaxis()->GetXmax();
		std::vector<double> binEdgesY(nBinsY+1);
		if(useLogYBins && !setLogYOnPads){
    			std::cout << "Error: request to set logybins but not to use log-scale on the y axis is probably not intended" << std::endl;
    			std::_Exit(99);
		}
		if(useLogYBins){
			HelperFunctions::equidistLogBins(binEdgesY, ymin, ymax, useLogYBins);
			data->SetBins(nBinsX,&binEdgesX[0],nBinsY,&binEdgesY[0]);
		}
	}
}

void TopMassControlPlots::MyHistogram::ConfigurePlotRanges(double PlotRangeYRatioMin, double PlotRangeYRatioMax, double PlotRangeYMin, double PlotRangeYMax){
	plotRangeYMin = PlotRangeYMin;
	plotRangeYMax = PlotRangeYMax;
	plotRangeYRatioMin = PlotRangeYRatioMin;
	plotRangeYRatioMax = PlotRangeYRatioMax;
}

void TopMassControlPlots::MyHistogram::ConfigureMoreExtraOptions(bool plotStackNormalized){
  plotStackNorm = plotStackNormalized;
}

void TopMassControlPlots::MyHistogram::SetExtraLabel(std::string label){
  extraLabel = label;
}
