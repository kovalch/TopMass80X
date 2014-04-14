#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerJetMatch.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerJetMatch::AnalyzerJetMatch(const std::vector<TString>& selectionStepsNoCategories,
                                   const std::vector<TString>& stepsForCategories,
                                   const JetCategories* jetCategories):
AnalyzerBase("jetMatch_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up jet match analyzer\n";
    std::cout<<"=== Finishing setting up jet match analyzer\n\n";
}



void AnalyzerJetMatch::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "matchedBjetFromTop";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "matchedBjetFromTop;;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "bQuark-genJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "genJet-recoJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Matched");
    
    name = "unmatchedGenBjetFromTop";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "unmatchedGenBjetFromTop;;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Several hadrons");
    
    name = "unmatchedRecoBjetFromTop";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "unmatchedRecoBjetFromTop;;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    
    name = "matchedBjetFromHiggs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "matchedBjetFromHiggs;;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "bQuark-genJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "genJet-recoJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Matched");
    
    name = "unmatchedGenBjetFromHiggs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "unmatchedGenBjetFromHiggs;;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Higgs jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Several hadrons");
    
    name = "unmatchedRecoBjetFromHiggs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "unmatchedRecoBjetFromHiggs;;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Higgs jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    
    name = "matchedBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "matchedBjet;;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top-Higgs gen overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "Unique gen jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Top-Higgs reco overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Unique reco jets");
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    int bins_p            = 30;  int max_p            = 150;
    int bins_pt           = 50;  int max_pt           = 300;
    int bins_ratio_pt     = 30;  int max_ratio_pt     = 3;
    int bins_dif_pt       = 30;  int min_dif_pt       = -100;  int max_dif_pt         = 200;
    int bins_dif_ratio_pt = 120; float min_dif_ratio_pt =-1.2; float max_dif_ratio_pt = 1.2;

    int bins_eta          = 50;  float min_eta        = -2.6;  float max_eta          = 2.6;
    int bins_phi          = 40;  int min_phi          = -4;    int max_phi            = 4;

    int bins_MinDR_match  = 250;  int min_MinDR_match = 0;    int max_MinDR_match     = 5.;
    int bins_MinDR        = 260;  int min_MinDR        = 0;    int max_MinDR          = 2.;    

    // Gen-Reco matching (for all the jets no matter if they are Top or Higgs Jets)
    
    name = "MinDeltaR_RecoGen";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen-Reco Matching;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));

    name = "RecoPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));
    
    name = "GenPtMinusRecoPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));
    
    name = "GenPtVSRecoPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "GenEtaVSRecoEta_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "GenPhiVSRecoPhi_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "GenP_VS_RecoP_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));

    name = "GenPtVSGenPtMinusRecoPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "MinDeltaR_RecoGen_VS_PtGen_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Gen-Reco Matching;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 
    
    name = "RecoPtOverGenPtVsGenPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenEta_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenPhi_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsMinDeltaR_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));


     // WRONG Gen-Reco matching (for all the jets no matter if they are Top or Higgs Jets)
    
    name = "Mismatched_MinDeltaR_RecoGen";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen-Reco Matching;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));

    name = "Mismatched_RecoPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));
    
    name = "Mismatched_GenPtMinusRecoPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));
    
    name = "Mismatched_GenPtVSRecoPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "Mismatched_GenEtaVSRecoEta_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "Mismatched_GenPhiVSRecoPhi_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "Mismatched_GenP_VS_RecoP_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));

    name = "Mismatched_GenPtVSGenPtMinusRecoPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Gen-Reco Matching;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 
    
    name = "Mismatched_RecoPtOverGenPtVsGenPt_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenEta_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenPhi_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

    //==================================               Top jets              =========================================//
    
    name = "MinDeltaR_RecoGen_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jets;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));

    name = "RecoPtOverGenPt_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));
    
    name = "GenPtMinusRecoPtOverGenPt_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));
    
    name = "GenPtVSRecoPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "GenEtaVSRecoEta_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "GenPhiVSRecoPhi_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "GenP_VS_RecoP_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));

    name = "GenPtVSGenPtMinusRecoPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jets;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 
    
    name = "RecoPtOverGenPtVsGenPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenEta_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenPhi_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

    
    // WRONG Gen-Reco mathching
    
    name = "Mismatched_MinDeltaR_RecoGen_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jets;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));

    name = "Mismatched_RecoPtOverGenPt_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));
    
    name = "Mismatched_GenPtMinusRecoPtOverGenPt_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));
    
    name = "Mismatched_GenPtVSRecoPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "Mismatched_GenEtaVSRecoEta_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "Mismatched_GenPhiVSRecoPhi_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "Mismatched_GenP_VS_RecoP_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));

    name = "Mismatched_GenPtVSGenPtMinusRecoPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jets;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 
    
    name = "Mismatched_RecoPtOverGenPtVsGenPt_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenEta_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenPhi_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Top jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

   //==================================               Higgs Jets              =========================================//
    
    name = "MinDeltaR_RecoGen_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jets;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));
    
    name = "RecoPtOverGenPt_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));

    name = "GenPtMinusRecoPtOverGenPt_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

    name = "GenPtVSRecoPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "GenEtaVSRecoEta_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "GenPhiVSRecoPhi_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "GenP_VS_RecoP_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));
    
    name = "GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen}",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jets;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 

    name = "RecoPtOverGenPtVsGenPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenEta_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

    // WRONG Gen-Reco matching
    
    name = "Mismatched_MinDeltaR_RecoGen_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jets;Min_{#DeltaR_{gen,reco}};Entries",bins_MinDR_match,min_MinDR_match,max_MinDR_match));
    
    name = "Mismatched_RecoPtOverGenPt_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jet;p_{T,Reco}/p_{T,Gen};Entries",bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_GenPtMinusRecoPtOverGenPt_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jet;(p_{T,Gen}-p_{T,Reco})/p_{T,Gen};Entries",bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

    name = "Mismatched_GenPtVSRecoPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Reco} [GeV];p_{T,Gen} [GeV]",bins_pt,0,max_pt,bins_pt,0,max_pt));

    name = "Mismatched_GenEtaVSRecoEta_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#eta_{Reco};#eta_{Gen}",bins_eta,min_eta,max_eta,bins_eta,min_eta,max_eta));

    name = "Mismatched_GenPhiVSRecoPhi_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#phi_{Reco};#phi_{Gen}",bins_phi,min_phi,max_phi,bins_phi,min_phi,max_phi));

    name = "Mismatched_GenP_VS_RecoP_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;|p_{Reco}|;|p_{Gen}|",bins_p,0,max_p,bins_p,0,max_p));
    
    name = "Mismatched_GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Gen}-p_{T,Reco} [GeV];p_{T,Gen}",bins_dif_pt,min_dif_pt,max_dif_pt,bins_pt,0,max_pt));

    name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jets;p_{T,Gen} [GeV];Min_{#DeltaR_{gen,reco}}",bins_pt,0,max_pt,bins_MinDR,min_MinDR,max_MinDR)); 

    name = "Mismatched_RecoPtOverGenPtVsGenPt_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;p_{T,Gen} [GeV];p_{T,Reco}/p_{T,Gen}",bins_pt,0,max_pt,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenEta_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#eta_{Gen};p_{T,Reco}/p_{T,Gen}",bins_eta,min_eta,max_eta,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;#phi_{Gen};p_{T,Reco}/p_{T,Gen}",bins_phi,min_phi,max_phi,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;Min_{#DeltaR_{gen,reco}};p_{T,Reco}/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_ratio_pt,0,max_ratio_pt));

    name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Higgs jet;Min_{#DeltaR_{gen,reco}};(p_{T,Gen}-p_{T,Reco})/p_{T,Gen}",bins_MinDR_match,min_MinDR_match,max_MinDR_match,bins_dif_ratio_pt,min_dif_ratio_pt,max_dif_ratio_pt));

}



void AnalyzerJetMatch::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                  const KinRecoObjects&,
                                  const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                  const double& weight, const TString&,
                                  std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "matchedBjetFromTop";
    if(!genObjectIndices.uniqueGenTopMatching()) m_histogram[name]->Fill(0., weight);
    else if(!genObjectIndices.uniqueRecoTopMatching()) m_histogram[name]->Fill(1., weight);
    else m_histogram[name]->Fill(2., weight);

    name = "unmatchedGenBjetFromTop";
    if(genObjectIndices.genBjetFromTopIndex_==-2 || genObjectIndices.genAntiBjetFromTopIndex_==-2) m_histogram[name]->Fill(3., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==-1 && genObjectIndices.genAntiBjetFromTopIndex_==-1) m_histogram[name]->Fill(1., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==-1 || genObjectIndices.genAntiBjetFromTopIndex_==-1) m_histogram[name]->Fill(2., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==genObjectIndices.genAntiBjetFromTopIndex_) m_histogram[name]->Fill(0., weight);

    name = "unmatchedRecoBjetFromTop";
    if(genObjectIndices.uniqueGenTopMatching()){
        if(genObjectIndices.recoBjetFromTopIndex_<0 && genObjectIndices.recoAntiBjetFromTopIndex_<0) m_histogram[name]->Fill(1., weight);
        else if(genObjectIndices.recoBjetFromTopIndex_<0 || genObjectIndices.recoAntiBjetFromTopIndex_<0) m_histogram[name]->Fill(2., weight);
        else if(genObjectIndices.recoBjetFromTopIndex_==genObjectIndices.recoAntiBjetFromTopIndex_) m_histogram[name]->Fill(0., weight);
    }

    name = "matchedBjetFromHiggs";
    if(!genObjectIndices.uniqueGenHiggsMatching()) m_histogram[name]->Fill(0., weight);
    else if(!genObjectIndices.uniqueRecoHiggsMatching()) m_histogram[name]->Fill(1., weight);
    else m_histogram[name]->Fill(2., weight);

    name = "unmatchedGenBjetFromHiggs";
    if(genObjectIndices.genBjetFromHiggsIndex_==-2 || genObjectIndices.genAntiBjetFromHiggsIndex_==-2) m_histogram[name]->Fill(3., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==-1 && genObjectIndices.genAntiBjetFromHiggsIndex_==-1) m_histogram[name]->Fill(1., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==-1 || genObjectIndices.genAntiBjetFromHiggsIndex_==-1) m_histogram[name]->Fill(2., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==genObjectIndices.genAntiBjetFromHiggsIndex_) m_histogram[name]->Fill(0., weight);

    name = "unmatchedRecoBjetFromHiggs";
    if(genObjectIndices.uniqueGenHiggsMatching()){
        if(genObjectIndices.recoBjetFromHiggsIndex_<0 && genObjectIndices.recoAntiBjetFromHiggsIndex_<0) m_histogram[name]->Fill(1., weight);
        else if(genObjectIndices.recoBjetFromHiggsIndex_<0 || genObjectIndices.recoAntiBjetFromHiggsIndex_<0) m_histogram[name]->Fill(2., weight);
        else if(genObjectIndices.recoBjetFromHiggsIndex_==genObjectIndices.recoAntiBjetFromHiggsIndex_) m_histogram[name]->Fill(0., weight);
    }

    name = "matchedBjet";
    if(genObjectIndices.uniqueGenMatching()){
        m_histogram[name]->Fill(1., weight);
        if(genObjectIndices.uniqueRecoMatching()) m_histogram[name]->Fill(3., weight);
        else if(genObjectIndices.uniqueRecoTopMatching() && genObjectIndices.uniqueRecoHiggsMatching()) m_histogram[name]->Fill(2., weight);
    }
    else if(genObjectIndices.uniqueGenTopMatching() && genObjectIndices.uniqueGenHiggsMatching()) m_histogram[name]->Fill(0., weight);


      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
      
        
    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& genJets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
    std::vector<int> genJetsIndices = common::initialiseIndices(genJets);
    common::orderIndices(genJetsIndices, genJets, common::LVpt);


    std::vector<int> recoIndex;       recoIndex.clear();
    std::vector<int> topRecoIndex;    topRecoIndex.clear();
    std::vector<int> higgsRecoIndex;  higgsRecoIndex.clear();
    std::vector<int> genIndex;        genIndex.clear();
    std::vector<int> topGenIndex;     topGenIndex.clear();
    std::vector<int> higgsGenIndex;   higgsGenIndex.clear();
    std::vector<int> bHadGenJetIndex; bHadGenJetIndex.clear();

    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);

    double deltaR=999.; 
    double minR = 999.;              std::vector<double> minRRecoGen;      minRRecoGen.clear();
    double minRTopRecoGen_ = 999.;   std::vector<double> minRTopRecoGen;   minRTopRecoGen.clear();
    double minRHiggsRecoGen_ = 999.; std::vector<double> minRHiggsRecoGen; minRHiggsRecoGen.clear();

    bool isTopOrHiggsMatched = false; bool isTopMatched = false; bool isHiggsMatched = false; 
        
    //Loop over the GenJets and choose the jets that are also included in the Hadronic List (bHadJetIndex)
    if (genJetsIndices.size()>0){
        for(const int Jetindex : genJetsIndices){ 
            for (const int Hadindex : bHadJetIndex){
                if (Jetindex!=Hadindex) continue;
                bHadGenJetIndex.push_back(Hadindex); 
            }
        }
    }
        
    //============================     Find the Top and Higgs Jets  (Gen Level)  ================================================//
    std::vector<int> topGenJetsIndex;       topGenJetsIndex.clear();
    std::vector<int> higgsGenJetsIndex;     higgsGenJetsIndex.clear();
    std::vector<int> topHiggsGenJetsIndex;  topHiggsGenJetsIndex.clear();

    
    //Find the Top Jets
    if (genObjectIndices.genBjetFromTopIndex_>=0){
        topGenJetsIndex.push_back(genObjectIndices.genBjetFromTopIndex_);
        isTopMatched = true;
    }
    if (genObjectIndices.genAntiBjetFromTopIndex_>=0){
        topGenJetsIndex.push_back(genObjectIndices.genAntiBjetFromTopIndex_);
        isTopMatched = true;
    }
    
    //Find the Higgs Jets
    if (genObjectIndices.genBjetFromHiggsIndex_>=0){
        higgsGenJetsIndex.push_back(genObjectIndices.genBjetFromHiggsIndex_);
        isHiggsMatched = true;
    }
    if (genObjectIndices.genAntiBjetFromHiggsIndex_>=0){
        higgsGenJetsIndex.push_back(genObjectIndices.genAntiBjetFromHiggsIndex_);
        isHiggsMatched = true;
    }
    

    //Find Jets with Top and Higgs Jets
    for (size_t i=0;i!=topGenJetsIndex.size();i++){
        topHiggsGenJetsIndex.push_back(topGenJetsIndex.at(i));
    }
    for (size_t i=0;i!=higgsGenJetsIndex.size();i++){
        topHiggsGenJetsIndex.push_back(higgsGenJetsIndex.at(i));
    }
    
    if (isTopMatched || isHiggsMatched) isTopOrHiggsMatched = true;
    common::orderIndices(topGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
    common::orderIndices(higgsGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
    common::orderIndices(topHiggsGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
    

    //========================================    Gen-Reco matching    ================================================//
    //For each gen jet find the reco jet with the minimum DeltaR!!
    for (size_t i_gen=0;i_gen!=genJetsIndices.size();++i_gen){
        minR = 999.;
        bool isRecoToGenJetsMatched = false;
        int recoJetIndex(-1);
        int genJetIndex(-1);
        int recoTopJetIndex(-1);
        int genTopJetIndex(-1);
        int recoHiggsJetIndex(-1);
        int genHiggsJetIndex(-1);
       

        for(const int index : recoObjectIndices.jetIndices_){
            deltaR=999.;
            deltaR = ROOT::Math::VectorUtil::DeltaR(recoObjects.jets_->at(index),genJets.at(i_gen));
            
            //take only matched pairs
            if (deltaR>0.4) continue;
            if ( ((commonGenObjects.allGenJets_->at(i_gen).Pt()-recoObjects.jets_->at(index).Pt())/commonGenObjects.allGenJets_->at(i_gen).Pt())<-0.4) continue;
            if ( ((commonGenObjects.allGenJets_->at(i_gen).Pt()-recoObjects.jets_->at(index).Pt())/commonGenObjects.allGenJets_->at(i_gen).Pt())>0.6) continue;
            //if (deltaR<0.5) continue;  //test to check the events with the bad gen-reco matching
            if (deltaR==999.) std::cout<<"ERROR!!"<<std::endl;
             
            //choose smallest R
            if (deltaR>minR) continue;
            if (deltaR<minR){
                minR = deltaR;

                //grep the indices associated to the smallest R (No matter if the Gen Jets are Top or Higgs)
                genJetIndex = genJetsIndices.at(i_gen);
                recoJetIndex = index;
                isRecoToGenJetsMatched = true;
              
                //grep the indices associated to the smallest R (for the Top Gen Jets)
                if(topGenJetsIndex.size()>0){
                    for (size_t i_top=0;i_top!=topGenJetsIndex.size();i_top++){  //Check if the generated jet is a Top. If yes find the matching reco
                            if (topGenJetsIndex.at(i_top) == genJetsIndices.at(i_gen)){
                            genTopJetIndex = genJetsIndices.at(i_gen); 
                            recoTopJetIndex = index; 
                            minRTopRecoGen_ = minR; 
                            }  
                    }
                }
              
                //grep the indices associated to the smallest R (for the Higgs Gen Jets)
                if(higgsGenJetsIndex.size()>0){
                    for (size_t i_higgs=0;i_higgs!=higgsGenJetsIndex.size();i_higgs++){ //Check if the generated jet is a Higgs. If yes find the matching reco
                        if (higgsGenJetsIndex.at(i_higgs) == genJetsIndices.at(i_gen)){
                            genHiggsJetIndex = genJetsIndices.at(i_gen);
                            recoHiggsJetIndex = index;
                            minRHiggsRecoGen_ = minR; 
                        }  
                    }
                } 
            }
        }
         
        //avoid cases in which we have no smallest deltaR
        if (!isRecoToGenJetsMatched) continue;
        
        if (genJetIndex>=0 && recoJetIndex>=0){
                genIndex.push_back(genJetIndex);
                recoIndex.push_back(recoJetIndex);
                minRRecoGen.push_back(minR);
        }
        
        if (!isTopOrHiggsMatched) continue; 

        if (genTopJetIndex>=0 && recoTopJetIndex>=0){
                topGenIndex.push_back(genTopJetIndex);
                topRecoIndex.push_back(recoTopJetIndex);
                minRTopRecoGen.push_back(minRTopRecoGen_);
        }

        if (genHiggsJetIndex>=0 && recoHiggsJetIndex>=0){
                higgsGenIndex.push_back(genHiggsJetIndex);
                higgsRecoIndex.push_back(recoHiggsJetIndex);
                minRHiggsRecoGen.push_back(minRHiggsRecoGen_);
        }

    }  

    //============================  Find the events with the Mismatched Jets ==========================//

    //all jets (No matter if they are Top/Higgs/...)
    bool isInMismatchedJetsIndx = false;
    std::vector<bool> isInMismatchedJetsIndex; isInMismatchedJetsIndex.clear();
    if(genIndex.size()>1){
        for(size_t i=0;i!=genIndex.size();i++){
            isInMismatchedJetsIndx = MismatchedJets(genIndex, recoIndex, i);
            isInMismatchedJetsIndex.push_back(isInMismatchedJetsIndx);
        }
    }

    //Top Jets
    bool isInMismatchedTopJetsIndx = false;
    std::vector<bool> isInMismatchedTopJetsIndex; isInMismatchedTopJetsIndex.clear();
    if(topGenIndex.size()>1){
        for(size_t i=0;i!=topGenIndex.size();i++){
            isInMismatchedTopJetsIndx = MismatchedJets(topGenIndex, topRecoIndex, i);
            isInMismatchedTopJetsIndex.push_back(isInMismatchedTopJetsIndx);
        }
    }

    //Higgs Jets
    bool isInMismatchedHiggsJetsIndx = false;
    std::vector<bool> isInMismatchedHiggsJetsIndex; isInMismatchedHiggsJetsIndex.clear();
    if(higgsGenIndex.size()>1){
        for(size_t i=0;i!=higgsGenIndex.size();i++){
            isInMismatchedHiggsJetsIndx = MismatchedJets(higgsGenIndex, higgsRecoIndex, i);
            isInMismatchedHiggsJetsIndex.push_back(isInMismatchedHiggsJetsIndx);
        }
    }
 
    //=======================     Find if there are Overlapping hadrons   ==========================//
    bool isOverlapping = false;
    if(bHadGenJetIndex.size()>1){ //Count the overlapping hadrons only for the cases with at least 2 Jets
        bool nOverlapping = overlappingHadrons( bHadGenJetIndex);
        if (nOverlapping>0) isOverlapping=true;
        else if (nOverlapping==0) isOverlapping=false;
    }
    
      
    //===========================           Fill the histograms        ==============================//
    if(topGenIndex.size()>0 && genObjectIndices.uniqueGenTopMatching()){
        for(size_t iJet=0;iJet!=topGenIndex.size();iJet++){
            if (isInMismatchedTopJetsIndex.size()>0){
                if(isInMismatchedTopJetsIndex.at(iJet)) continue;
                double TopgenP  = commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).M();
                double ToprecoP = recoObjects.jets_->at(topRecoIndex.at(iJet)).M();
          
                name = "GenP_VS_RecoP_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(ToprecoP,TopgenP,weight);
          
                name = "MinDeltaR_RecoGen_Topjet";
                m_histogram[name]->Fill(minRTopRecoGen.at(iJet),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(), minRTopRecoGen.at(iJet),weight);

                name = "RecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtVSRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);
      
                name = "GenEtaVSRecoEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Eta(),weight);

                name = "GenPhiVSRecoPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Eta(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Phi(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(iJet), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);
            }
        }
    }

    if(higgsGenIndex.size()>0 && genObjectIndices.uniqueGenHiggsMatching()){
        for(size_t iJet=0;iJet!=higgsGenIndex.size();iJet++){
            if (isInMismatchedHiggsJetsIndex.size()>0){
                if(isInMismatchedHiggsJetsIndex.at(iJet)) continue;
                double HiggsgenP  = commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).M();
                double HiggsrecoP = recoObjects.jets_->at(higgsRecoIndex.at(iJet)).M();          
                
                name = "GenP_VS_RecoP_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(HiggsrecoP,HiggsgenP,weight);
                
                name = "MinDeltaR_RecoGen_Higgsjet";
                m_histogram[name]->Fill(minRHiggsRecoGen.at(iJet),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(), minRHiggsRecoGen.at(iJet),weight);

                name = "RecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtVSRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);
          
                name = "GenEtaVSRecoEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Eta(),weight);

                name = "GenPhiVSRecoPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Eta(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Phi(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(iJet), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);
            }
        }
    }


    if(genIndex.size()>0 && !isOverlapping){
        for(size_t iJet=0;iJet!=genIndex.size();iJet++){
            if (isInMismatchedJetsIndex.size()>0){
                if(isInMismatchedJetsIndex.at(iJet)) continue;
                double genP  = commonGenObjects.allGenJets_->at(genIndex.at(iJet)).M();
                double recoP = recoObjects.jets_->at(recoIndex.at(iJet)).M();          
                
                name = "GenP_VS_RecoP_2D";
                ((TH2D*)m_histogram[name])->Fill(recoP,genP,weight);
                
                name = "MinDeltaR_RecoGen";
                m_histogram[name]->Fill(minRRecoGen.at(iJet),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(), minRRecoGen.at(iJet),weight);

                name = "RecoPtOverGenPt";
                m_histogram[name]->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "GenPtVSRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);
          
                name = "GenEtaVSRecoEta_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Eta(),weight);

                name = "GenPhiVSRecoPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Eta(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Phi(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(iJet), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);
            }
        }
    }

    //Plots for the WRONG Matching
    if(genIndex.size()>0 && !isOverlapping){
        for(size_t iJet=0;iJet!=genIndex.size();iJet++){
            if (isInMismatchedJetsIndex.size()>0){
                if(!isInMismatchedJetsIndex.at(iJet)) continue;
                double genP  = commonGenObjects.allGenJets_->at(genIndex.at(iJet)).M();
                double recoP = recoObjects.jets_->at(recoIndex.at(iJet)).M();  
                
                name = "Mismatched_GenP_VS_RecoP_2D";
                ((TH2D*)m_histogram[name])->Fill(recoP,genP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen";
                m_histogram[name]->Fill(minRRecoGen.at(iJet),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(), minRRecoGen.at(iJet),weight);

                name = "Mismatched_RecoPtOverGenPt";
                m_histogram[name]->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Eta(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Phi(), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(iJet), recoObjects.jets_->at(recoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt()-recoObjects.jets_->at(recoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(iJet)).Pt(),weight);
            }
        }
    } 
    
    if(topGenIndex.size()>0  && genObjectIndices.uniqueGenTopMatching()){
        for(size_t iJet=0;iJet!=topGenIndex.size();iJet++){
            if (isInMismatchedTopJetsIndex.size()>0){
                if(!isInMismatchedTopJetsIndex.at(iJet)) continue;
                double TopgenP  = commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).M();
                double ToprecoP = recoObjects.jets_->at(topRecoIndex.at(iJet)).M();
                
                name = "Mismatched_GenP_VS_RecoP_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(ToprecoP,TopgenP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen_Topjet";
                m_histogram[name]->Fill(minRTopRecoGen.at(iJet),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(), minRTopRecoGen.at(iJet),weight);

                name = "Mismatched_RecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(topRecoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Eta(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Phi(), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(iJet), recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(topRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(topGenIndex.at(iJet)).Pt(),weight);
            }
        }
    }

    if(higgsGenIndex.size()>0 && genObjectIndices.uniqueGenHiggsMatching()){
        for(size_t iJet=0;iJet!=higgsGenIndex.size();iJet++){
            if (isInMismatchedHiggsJetsIndex.size()>0){
                if(!isInMismatchedHiggsJetsIndex.at(iJet)) continue;
                double HiggsgenP  = commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).M();
                double HiggsrecoP = recoObjects.jets_->at(higgsRecoIndex.at(iJet)).M();  
                
                name = "Mismatched_GenP_VS_RecoP_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(HiggsrecoP,HiggsgenP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen_Higgsjet";
                m_histogram[name]->Fill(minRHiggsRecoGen.at(iJet),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(), minRHiggsRecoGen.at(iJet),weight);

                name = "Mismatched_RecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Eta(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Phi(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt(), commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Eta(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Phi(), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(iJet), recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt()/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(iJet), (commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt()-recoObjects.jets_->at(higgsRecoIndex.at(iJet)).Pt())/commonGenObjects.allGenJets_->at(higgsGenIndex.at(iJet)).Pt(),weight);
            }
        }
    }

    
}


int AnalyzerJetMatch::overlappingHadrons( const std::vector<int>& genIndx)
{
    
    bool isInOverlappingJetsIndex = false; bool isInOverlappingGenJetsIndex = false;
    int n_Overlapping=0;
    std::vector<int> OverlappingIndex; OverlappingIndex.clear();
    std::vector<int> OverlappingGenIndex; OverlappingGenIndex.clear();
    
    //Count how many (in total) overlapping hadrons there are in each event
    for(size_t i=0;i!=genIndx.size();i++){
        for(size_t j=i+1;j!=genIndx.size();j++){
            if(genIndx.at(i)==genIndx.at(j)){
                if (OverlappingIndex.size()<1){
                    OverlappingGenIndex.push_back(genIndx.at(j));
                    OverlappingIndex.push_back(j);
                }  
                else if (OverlappingIndex.size()>0){
                    for (size_t k=0; k!=OverlappingIndex.size(); k++){ 
                        if (static_cast<int>(j)==OverlappingIndex.at(k)){
                            isInOverlappingJetsIndex = true;
                            break;
                        }
                        else if (static_cast<int>(j)!=OverlappingIndex.at(k)) isInOverlappingJetsIndex = false; 
                    }
                    for (size_t g=0; g!=OverlappingGenIndex.size(); g++){
                        if (genIndx.at(j)==OverlappingGenIndex.at(g)){ 
                            isInOverlappingGenJetsIndex = true;
                            break;
                        }
                        else if (genIndx.at(j)!=OverlappingGenIndex.at(g)) isInOverlappingGenJetsIndex = false;
                    }
                    
                    
                    if (!isInOverlappingJetsIndex) OverlappingIndex.push_back(j); 
                    if (!isInOverlappingGenJetsIndex)  OverlappingGenIndex.push_back(genIndx.at(j));
                }
                
                if (OverlappingIndex.size()>0){
                    for(size_t n=0; n!=OverlappingIndex.size();n++){
                        if(static_cast<int>(i)!=OverlappingIndex.at(n)){
                            isInOverlappingJetsIndex = false;
                        }
                        else if(static_cast<int>(i)==OverlappingIndex.at(n)){
                            isInOverlappingJetsIndex = true;
                            break;
                        }  
                    }
                }            
                
                if (!isInOverlappingJetsIndex){ 
                    n_Overlapping=OverlappingIndex.size()+OverlappingGenIndex.size();
                }
            }
            
        }
    }

    return n_Overlapping;
    
}



bool AnalyzerJetMatch::MismatchedJets( const std::vector<int>& genIndx, const std::vector<int>& recoIndx, const int& i)
{
    
    bool inMismatchedJetIndex = false;
    
    std::vector<int> MismatchedIndex; MismatchedIndex.clear();
    
    for(size_t j=i+1;j!=genIndx.size();j++){       
        if(recoIndx.at(i)==recoIndx.at(j) && genIndx.at(i)!=genIndx.at(j)){
            
            if(MismatchedIndex.size()<1){
                inMismatchedJetIndex = true;
                MismatchedIndex.push_back(i);
                MismatchedIndex.push_back(j);
            }
            
            else if (MismatchedIndex.size()>0){
                for (size_t k=0; k!=MismatchedIndex.size(); k++){
                    if (static_cast<int>(i)!=MismatchedIndex.at(k)) inMismatchedJetIndex = false;
                    else if (static_cast<int>(i)==MismatchedIndex.at(k)){
                        inMismatchedJetIndex = true;
                        break;
                    }
                }
                
                if (!inMismatchedJetIndex) MismatchedIndex.push_back(i);   
                
                for (size_t k=0; k!=MismatchedIndex.size(); k++){
                    if (static_cast<int>(j)!=MismatchedIndex.at(k)) inMismatchedJetIndex = false; 
                    else if (static_cast<int>(j)==MismatchedIndex.at(k)){
                        inMismatchedJetIndex = true;
                        break;
                    }
                }
                
                if (!inMismatchedJetIndex) MismatchedIndex.push_back(j);   
            }
            
        }   
    }
    if (MismatchedIndex.size()>0){
        for(size_t n=0; n!=MismatchedIndex.size();n++){
            if(static_cast<int>(i)!=MismatchedIndex.at(n)){
                inMismatchedJetIndex = false;
            }
            else if(static_cast<int>(i)==MismatchedIndex.at(n)){
                inMismatchedJetIndex = true;
                break;
            }  
        }
    }
    return inMismatchedJetIndex;
}

