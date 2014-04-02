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

    int bins_DEta         = 50;  int min_DEta         = -5;    int max_DEta           = 5;
    int bins_MinDEta         = 40;  int min_MinDEta   = 0;   int max_MinDEta          = 4;

    int bins_MinDPhi      = 32;  int min_MinDPhi      = 0;  int max_MinDPhi           = 3;
    int bins_DR           = 50;  int min_DR           = 0;    int max_DR              = 5;
    int bins_MinDR_match  = 250;  int min_MinDR_match = 0;    int max_MinDR_match     = 5.;
  
    int bins_MinDR        = 260;  int min_MinDR        = 0;    int max_MinDR          = 2.;    

    name = "NOverlappingHadrons_Top";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Overlapping Jets Multiplicity;# overlapping hadrons;Entries",10,0,10));
    
    name = "HasOverlappingHadrons_Top";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Has Overlapping Hadrons; overlapping hadrons;Entries",2,0,2));
     
    name = "GenJetsID_VS_NOverlappingHadrons_Top_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"# GenJets Vs # Overlapping Hadrons;# overlapping hadrons ;# Gen jets",10,0,10,5,0,5));
   

    name = "NOverlappingHadrons_Higgs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Overlapping Jets Multiplicity;# overlapping hadrons;Entries",10,0,10));
    
    name = "HasOverlappingHadrons_Higgs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Has Overlapping Hadrons; overlapping hadrons;Entries",2,0,2));
     
    name = "GenJetsID_VS_NOverlappingHadrons_Higgs_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"# GenJets Vs # Overlapping Hadrons;# overlapping hadrons ;# Gen jets",10,0,10,5,0,5));

    name = "NOverlappingHadrons_TopHiggs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Overlapping Jets Multiplicity;# overlapping hadrons;Entries",10,0,10));
    
    name = "HasOverlappingHadrons_TopHiggs";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Has Overlapping Hadrons; overlapping hadrons;Entries",2,0,2));
     
    name = "GenJetsID_VS_NOverlappingHadrons_TopHiggs_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"# GenJets Vs # Overlapping Hadrons;# overlapping hadrons ;# Gen jets",10,0,10,5,0,5));

     
    name = "NOverlappingHadrons";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Overlapping Jets Multiplicity;# overlapping hadrons;Entries",10,0,10));
    
    name = "HasOverlappingHadrons";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Has Overlapping Hadrons; overlapping hadrons;Entries",2,0,2));
     
    name = "GenJetsID_VS_NOverlappingHadrons_2D";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"# GenJets Vs # Overlapping Hadrons;# overlapping hadrons ;# Gen jets",10,0,10,5,0,5));
  
     
    name = "MinDeltaR";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs Jets;Min_{#DeltaR_{jj}};Entries",bins_MinDR,min_MinDR,max_MinDR));

    name = "MinDeltaEta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#eta_{jj}};Entries",bins_MinDEta,min_MinDEta,max_MinDEta));
    
    name = "MinDeltaPhi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#phi_{jj}};Entries",bins_MinDPhi,min_MinDPhi,max_MinDPhi));
    
    
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
    name = "Topleading_TopHiggsJets_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet in Top+Higgs Jets;p_{T} [GeV];Entries",bins_pt,0,max_pt));
    
    name = "leading_Topjet_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet;p_{T} [GeV];Entries",bins_pt,0,max_pt));

    name = "Topleading_TopHiggsJets_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet in Top+Higgs Jets;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "leading_Topjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "Topsubleading_TopHiggsJets_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet in Top+Higgs Jets;p_{T} [GeV];Entries",bins_pt,0,max_pt));

    name = "subleading_Topjet_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet;p_{T} [GeV];Entries",bins_pt,0,max_pt));
    
    name = "Topsubleading_TopHiggsJets_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet in Top+Higgs Jets;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "subleading_Topjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet;#eta;Entries",bins_eta,min_eta,max_eta));


    name = "DeltaEta_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jets;#Delta#eta_{jj};Entries",bins_DEta,min_DEta,max_DEta));
    
    name = "DeltaPhi_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jets;#Delta#phi_{jj};Entries",bins_phi,min_phi,max_phi));

    name = "DeltaR_Topjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top jets;#DeltaR_{jj};Entries",bins_DR,min_DR,max_DR));

        
    
    // Gen-Reco mathching
    
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
    name = "Higgsleading_TopHiggsJets_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet in Top+Higgs Jets;p_{T} [GeV];Entries",bins_pt,0,max_pt));

    name = "leading_Higgsjet_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet;p_{T} [GeV];Entries",bins_pt,0,max_pt));
    
    name = "Higgsleading_TopHiggsJets_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet in Top+Higgs Jets;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "leading_Higgsjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "Higgssubleading_TopHiggsJets_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subeading Higgs jet in Top+Higgs Jets;p_{T} [GeV];Entries",bins_pt,0,max_pt));

    name = "subleading_Higgsjet_Pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Higgs jet;p_{T} [GeV];Entries",bins_pt,0,max_pt));
    
    name = "Higgssubleading_TopHiggsJets_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Higgs jet in Top+Higgs Jets;#eta;Entries",bins_eta,min_eta,max_eta));

    name = "subleading_Higgsjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Higgs jet;#eta;Entries",bins_eta,min_eta,max_eta));

   
    name = "DeltaEta_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jets;#Delta#eta_{jj};Entries",bins_DEta,min_DEta,max_DEta));
    
    name = "DeltaPhi_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jets;#Delta#phi_{jj};Entries",bins_phi,min_phi,max_phi));

    name = "DeltaR_Higgsjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Higgs jets;#DeltaR_{jj};Entries",bins_DR,min_DR,max_DR));
    
    
     // Gen-Reco matching
    
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
    const VLV& testcommonjets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
    std::vector<int> testcommonJetIndices = common::initialiseIndices(testcommonjets);
    common::orderIndices(testcommonJetIndices, testcommonjets, common::LVpt);


    std::vector<int> recoIndex;       recoIndex.clear();
    std::vector<int> ToprecoIndex;    ToprecoIndex.clear();
    std::vector<int> HiggsrecoIndex;  HiggsrecoIndex.clear();
    std::vector<int> genIndex;        genIndex.clear();
    std::vector<int> TopgenIndex;     TopgenIndex.clear();
    std::vector<int> HiggsgenIndex;   HiggsgenIndex.clear();
    std::vector<int> flavour; std::vector<int> testbHadcommonJetIndex; 

    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);

    double deltaR=999.; double deltaEta=999.; double deltaPhi=999.;
                        double deltaEta_Abs=999.; double deltaPhi_Abs=999.9; 
    double minR = 999.; double minEta = 999.; double minPhi = 999.;
                                     std::vector<double> minRRecoGen; minRRecoGen.clear();
    double minRTopRecoGen_ = 999.;   std::vector<double> minRTopRecoGen;   minRTopRecoGen.clear();
    double minRHiggsRecoGen_ = 999.; std::vector<double> minRHiggsRecoGen; minRHiggsRecoGen.clear();

    bool isTopOrHiggsMatched = false; bool isTopMatched = false; bool isHiggsMatched = false;
    bool isOverlapping_Top = false; bool isOverlapping_Higgs = false; bool isOverlapping_TopHiggs = false; bool isOverlapping = false;
        
    //We loop over our jet list (testcommonjet) and choose the jets that are also included in the Hadronic List (bHadJetIndex). 

    if (testcommonJetIndices.size()>0){
        for(const int Jetindex : testcommonJetIndices){ 
            for (const int Hadindex : bHadJetIndex){
                if (Jetindex!=Hadindex) continue;
                testbHadcommonJetIndex.push_back(Hadindex); 
            }//for i_had
        }//for Jetindex
    }//if (testcommonJetIndices.size()>0)
        
    //============================     Find the Top and Higgs Jets  (Gen Level)  ================================================//
    std::vector<int> MyTopJetsIndex;       MyTopJetsIndex.clear();
    std::vector<int> MyHiggsJetsIndex;     MyHiggsJetsIndex.clear();
    std::vector<int> MyTopHiggsJetsIndex;  MyTopHiggsJetsIndex.clear();


    //Find the Top Jets
    if (genObjectIndices.genBjetFromTopIndex_>=0){
        MyTopJetsIndex.push_back(genObjectIndices.genBjetFromTopIndex_);
        isTopMatched = true;
    }
    if (genObjectIndices.genAntiBjetFromTopIndex_>=0){
        MyTopJetsIndex.push_back(genObjectIndices.genAntiBjetFromTopIndex_);
        isTopMatched = true;
    }
    
    //Find the Higgs Jets
    if (genObjectIndices.genBjetFromHiggsIndex_>=0){
        MyHiggsJetsIndex.push_back(genObjectIndices.genBjetFromHiggsIndex_);
        isHiggsMatched = true;
    }
    if (genObjectIndices.genAntiBjetFromHiggsIndex_>=0){
        MyHiggsJetsIndex.push_back(genObjectIndices.genAntiBjetFromHiggsIndex_);
        isHiggsMatched = true;
    }
 

    //Find Jets with Top and Higgs
    for (size_t i=0;i!=MyTopJetsIndex.size();i++){
        MyTopHiggsJetsIndex.push_back(MyTopJetsIndex.at(i));
    }
    for (size_t i=0;i!=MyHiggsJetsIndex.size();i++){
        MyTopHiggsJetsIndex.push_back(MyHiggsJetsIndex.at(i));
    }

    if (isTopMatched || isHiggsMatched) isTopOrHiggsMatched = true;
    common::orderIndices(MyTopJetsIndex, testcommonjets, common::LVpt);  //pt ordered jets 
    common::orderIndices(MyHiggsJetsIndex, testcommonjets, common::LVpt);  //pt ordered jets 
    common::orderIndices(MyTopHiggsJetsIndex, testcommonjets, common::LVpt);  //pt ordered jets 

    //============================     Find the Top and Higgs Jets  (Reco Level)  ================================================//

    //========================================    Gen-Reco matching    ================================================//
    //For each gen jet we find the reco jet with the minimum DeltaR!!
    for (size_t i_gen=0;i_gen!=testcommonJetIndices.size();++i_gen){
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
            deltaR = ROOT::Math::VectorUtil::DeltaR(recoObjects.jets_->at(index),testcommonjets.at(i_gen));
            
            //only take matched pairs
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
                genJetIndex = testcommonJetIndices.at(i_gen);
                recoJetIndex = index;
                isRecoToGenJetsMatched = true;
              
                //grep the indices associated to the smallest R (for the Top Gen Jets)
                if(MyTopJetsIndex.size()>0){
                    for (size_t i_top=0;i_top!=MyTopJetsIndex.size();i_top++){  //here I check if the generated jet is a Top. If yes I also find the matching reco
                            if (MyTopJetsIndex.at(i_top) == testcommonJetIndices.at(i_gen)){
                            genTopJetIndex = testcommonJetIndices.at(i_gen); 
                            recoTopJetIndex = index; 
                            minRTopRecoGen_ = minR; 
                            }  
                    }// for i_top
                }
              
                //grep the indices associated to the smallest R (for the Higgs Gen Jets)
                if(MyHiggsJetsIndex.size()>0){
                    for (size_t i_higgs=0;i_higgs!=MyHiggsJetsIndex.size();i_higgs++){ //here I check if the generated jet is a Higgs. If yes I also find the matching reco
                        if (MyHiggsJetsIndex.at(i_higgs) == testcommonJetIndices.at(i_gen)){
                            genHiggsJetIndex = testcommonJetIndices.at(i_gen);
                            recoHiggsJetIndex = index;
                            minRHiggsRecoGen_ = minR; 
                        }  
                    }// for i_higgs 
                } 
            }//if (deltaR<minR)
        }//index (reco)
         
        //avoid cases in which we have no smallest deltaR
        if (!isRecoToGenJetsMatched) continue;
        
        if (genJetIndex>=0 && recoJetIndex>=0){
                genIndex.push_back(genJetIndex);
                recoIndex.push_back(recoJetIndex);
                minRRecoGen.push_back(minR);
        }
        
        if (!isTopOrHiggsMatched) continue; 

        if (genTopJetIndex>=0 && recoTopJetIndex>=0){
                TopgenIndex.push_back(genTopJetIndex);
                ToprecoIndex.push_back(recoTopJetIndex);
                minRTopRecoGen.push_back(minRTopRecoGen_);
        }

        if (genHiggsJetIndex>=0 && recoHiggsJetIndex>=0){
                HiggsgenIndex.push_back(genHiggsJetIndex);
                HiggsrecoIndex.push_back(recoHiggsJetIndex);
                minRHiggsRecoGen.push_back(minRHiggsRecoGen_);
        }

    }//i_jen  

    bool TopHasLeadingAndSubleading = true;
    bool HiggsHasLeadingAndSubleading = true;

    bool isInOverlappingTopJetsIndex = false;
    bool isInOverlappingTopGenJetsIndex = false;

    bool isInOverlappingHiggsJetsIndex = false;
    bool isInOverlappingHiggsGenJetsIndex = false;

    bool isInOverlappingTopHiggsJetsIndex = false;
    bool isInOverlappingTopHiggsGenJetsIndex = false;

    bool isInOverlappingJetsIndex = false;
    bool isInOverlappingGenJetsIndex = false;

    //Top Jets
    int LeadingTopJetIndx =-1;
    int SubLeadingTopJetIndx =-1;

    //Higgs Jets
    int LeadingHiggsJetIndx =-1;
    int SubLeadingHiggsJetIndx =-1;

    //=======================================  Find the events with the Mismatched Jets ===================================//
    //====  all jets (No matter if they are Top/Higgs/...)
    bool isInMismatchedJetsIndex_ = false;
    std::vector<bool> isInMismatchedJetsIndex; isInMismatchedJetsIndex.clear();

    std::vector<int> MismatchedIndex; MismatchedIndex.clear();
    std::vector<int> MismatchedGenIndex; MismatchedGenIndex.clear();
    std::vector<int> MismatchedRecoIndex; MismatchedRecoIndex.clear();
    std::vector<int> MismatchedDRIndex; MismatchedDRIndex.clear();
    
    if(genIndex.size()>1){
        for(size_t i=0;i!=genIndex.size();i++){
            for(size_t j=i+1;j!=genIndex.size();j++){       
                if(recoIndex.at(i)==recoIndex.at(j) && genIndex.at(i)!=genIndex.at(j)){
  
                    if(MismatchedIndex.size()<1){
                        isInMismatchedJetsIndex_ = true;
                        MismatchedIndex.push_back(i);
                        MismatchedIndex.push_back(j);
                    }
             
                    else if (MismatchedIndex.size()>0){
                        for (size_t k=0; k!=MismatchedIndex.size(); k++){
                            //if (i!=MismatchedIndex.at(k)) isInMismatchedJetsIndex_ = false;
                            if (static_cast<int>(i)!=MismatchedIndex.at(k)) isInMismatchedJetsIndex_ = false;
                            //else if (i==MismatchedIndex.at(k)){
                            else if (static_cast<int>(i)==MismatchedIndex.at(k)){
                                isInMismatchedJetsIndex_ = true;
                                break;
                            }
                        }// for k
                
                        if (!isInMismatchedJetsIndex_) MismatchedIndex.push_back(i);   
            
                        for (size_t k=0; k!=MismatchedIndex.size(); k++){
                            //if (j!=MismatchedIndex.at(k)) isInMismatchedJetsIndex_ = false; 
                            if (static_cast<int>(j)!=MismatchedIndex.at(k)) isInMismatchedJetsIndex_ = false; 
                            //else if (j==MismatchedIndex.at(k)){
                            else if (static_cast<int>(j)==MismatchedIndex.at(k)){
                                isInMismatchedJetsIndex_ = true;
                                break;
                            }
                        }// for k
                
                        if (!isInMismatchedJetsIndex_) MismatchedIndex.push_back(j);   
                    }//else if (MismatchedIndex.size()>0) 
            
                }//if(recoIndex.at(i)==recoIndex.at(j) && genIndex.at(i)!=genIndex.at(j))
                  
              
            }//for j
            if (MismatchedIndex.size()>0){
                for(size_t n=0; n!=MismatchedIndex.size();n++){
                    //if(i!=MismatchedIndex.at(n)){
                    if(static_cast<int>(i)!=MismatchedIndex.at(n)){
                        isInMismatchedJetsIndex_ = false;
                    }//if(i!=MismatchedIndex.at(n))
                    //else if(i==MismatchedIndex.at(n)){
                    else if(static_cast<int>(i)==MismatchedIndex.at(n)){
                        isInMismatchedJetsIndex_ = true;
                        break;
                    }//else if(i==MismatchedIndex.at(n))   
                }//for n 
            }//if (MismatchedIndex.size()>0)
            isInMismatchedJetsIndex.push_back(isInMismatchedJetsIndex_);        
        }//for i
    }//if(genIndex.size()>1)


    //===== Top Jets  ====//
    bool isInMismatchedTopJetsIndex_ = false;
    std::vector<bool> isInMismatchedTopJetsIndex; isInMismatchedTopJetsIndex.clear();

    std::vector<int> TopMismatchedIndex; TopMismatchedIndex.clear();
    std::vector<int> TopMismatchedGenIndex; TopMismatchedGenIndex.clear();
    std::vector<int> TopMismatchedRecoIndex; TopMismatchedRecoIndex.clear();
    std::vector<int> TopMismatchedDRIndex; TopMismatchedDRIndex.clear();
    //Correct it!!! And test it!!!
    if(TopgenIndex.size()>1){
        for(size_t i=0;i!=TopgenIndex.size();i++){
            for(size_t j=i+1;j!=TopgenIndex.size();j++){
                if(ToprecoIndex.at(i)==ToprecoIndex.at(j) && TopgenIndex.at(i)!=TopgenIndex.at(j)){
                    if(TopMismatchedIndex.size()<1){
                        isInMismatchedTopJetsIndex_ = true;
                        TopMismatchedIndex.push_back(i);          
                        TopMismatchedIndex.push_back(j);
                    }
             
                    else if (TopMismatchedIndex.size()>0){
                        for (size_t k=0; k!=TopMismatchedIndex.size(); k++){ 
                            //if (i!=TopMismatchedIndex.at(k)) isInMismatchedTopJetsIndex_ = false; 
                            if (static_cast<int>(i)!=TopMismatchedIndex.at(k)) isInMismatchedTopJetsIndex_ = false; 
                            //else if (i==TopMismatchedIndex.at(k)){
                            else if (static_cast<int>(i)==TopMismatchedIndex.at(k)){
                                isInMismatchedTopJetsIndex_ = true;
                                break;
                            }
                        }// for k
                
                        if (!isInMismatchedTopJetsIndex_) TopMismatchedIndex.push_back(i);   
            
                        for (size_t k=0; k!=TopMismatchedIndex.size(); k++){ 
                            //if (j!=TopMismatchedIndex.at(k)) isInMismatchedTopJetsIndex_ = false; 
                            //else if (j==TopMismatchedIndex.at(k)){
                            if (static_cast<int>(j)!=TopMismatchedIndex.at(k)) isInMismatchedTopJetsIndex_ = false; 
                            else if (static_cast<int>(j)==TopMismatchedIndex.at(k)){
                                isInMismatchedTopJetsIndex_ = true;
                                break;
                            }
                        }// for k
                
                        if (!isInMismatchedTopJetsIndex_) TopMismatchedIndex.push_back(j);   
                    }//else if (TopMismatchedIndex.size()>0) 
                }//if(ToprecoIndex.at(i)==ToprecoIndex.at(j) && TopgenIndex.at(i)!=TopgenIndex.at(j)) 
            }//for j
            if (TopMismatchedIndex.size()>0){
                for(size_t n=0; n!=TopMismatchedIndex.size();n++){
                    //if(i!=TopMismatchedIndex.at(n)){
                    if(static_cast<int>(i)!=TopMismatchedIndex.at(n)){
                        isInMismatchedTopJetsIndex_ = false;
                    }//if(i!=TopMismatchedIndex.at(n))
                    //else if(i==TopMismatchedIndex.at(n)){
                    else if(static_cast<int>(i)==TopMismatchedIndex.at(n)){
                        isInMismatchedTopJetsIndex_ = true;
                        break;
                    }//else if(i==TopMismatchedIndex.at(n))   
                }//for n 
            }//if (TopMismatchedIndex.size()>0)
            isInMismatchedTopJetsIndex.push_back(isInMismatchedTopJetsIndex_);        
        }//for i
    }//if


    //===== Higgs Jets  ====//
    bool isInMismatchedHiggsJetsIndex_ = false;
    std::vector<bool> isInMismatchedHiggsJetsIndex; isInMismatchedHiggsJetsIndex.clear();

    std::vector<int> HiggsMismatchedIndex; HiggsMismatchedIndex.clear();
    std::vector<int> HiggsMismatchedGenIndex; HiggsMismatchedGenIndex.clear();
    std::vector<int> HiggsMismatchedRecoIndex; HiggsMismatchedRecoIndex.clear();
    std::vector<int> HiggsMismatchedDRIndex; HiggsMismatchedDRIndex.clear();
    //Correct it!!! And test it!!!
    if(HiggsgenIndex.size()>1){
        for(size_t i=0;i!=HiggsgenIndex.size();i++){
            for(size_t j=i+1;j!=HiggsgenIndex.size();j++){
                    if(HiggsrecoIndex.at(i)==HiggsrecoIndex.at(j) && HiggsgenIndex.at(i)!=HiggsgenIndex.at(j)){
                        if(HiggsMismatchedIndex.size()<1){
                            isInMismatchedHiggsJetsIndex_ = true;
                            HiggsMismatchedIndex.push_back(i);  
                            HiggsMismatchedIndex.push_back(j);
                        }
                         
                        else if (HiggsMismatchedIndex.size()>0){
                            for (size_t k=0; k!=HiggsMismatchedIndex.size(); k++){ 
                                //if (i!=HiggsMismatchedIndex.at(k)) isInMismatchedHiggsJetsIndex_ = false; 
                                //else if (i==HiggsMismatchedIndex.at(k)){
                                if (static_cast<int>(i)!=HiggsMismatchedIndex.at(k)) isInMismatchedHiggsJetsIndex_ = false; 
                                else if (static_cast<int>(i)==HiggsMismatchedIndex.at(k)){
                                    isInMismatchedHiggsJetsIndex_ = true;
                                    break;
                                }
                            }// for k
                            
                            if (!isInMismatchedHiggsJetsIndex_) HiggsMismatchedIndex.push_back(i);   
                        
                            for (size_t k=0; k!=HiggsMismatchedIndex.size(); k++){
                                //if (j!=HiggsMismatchedIndex.at(k)) isInMismatchedHiggsJetsIndex_ = false; 
                                //else if (j==HiggsMismatchedIndex.at(k)){
                                if (static_cast<int>(j)!=HiggsMismatchedIndex.at(k)) isInMismatchedHiggsJetsIndex_ = false; 
                                else if (static_cast<int>(j)==HiggsMismatchedIndex.at(k)){
                                    isInMismatchedHiggsJetsIndex_ = true;
                                    break;
                                }
                            }// for k
                            
                            if (!isInMismatchedHiggsJetsIndex_) HiggsMismatchedIndex.push_back(j);   
                        }//else if (HiggsMismatchedIndex.size()>0)
                    }//if(HiggsrecoIndex.at(i)==HiggsrecoIndex.at(j) && HiggsgenIndex.at(i)!=HiggsgenIndex.at(j))
                              
                      
            }//for j
            if (HiggsMismatchedIndex.size()>0){
                for(size_t n=0; n!=HiggsMismatchedIndex.size();n++){
                    //if(i!=HiggsMismatchedIndex.at(n)){
                    if(static_cast<int>(i)!=HiggsMismatchedIndex.at(n)){
                            isInMismatchedHiggsJetsIndex_ = false;
                    }//if(i!=HiggsMismatchedIndex.at(n))
                    //else if(i==HiggsMismatchedIndex.at(n)){
                    else if(static_cast<int>(i)==HiggsMismatchedIndex.at(n)){
                            isInMismatchedHiggsJetsIndex_ = true;
                            break;
                    }//else if(i==HiggsMismatchedIndex.at(n))   
                }//for n 
            }//if (HiggsMismatchedIndex.size()>0)
            isInMismatchedHiggsJetsIndex.push_back(isInMismatchedHiggsJetsIndex_);
        }//for i
    }//if
   
    //==================================================================================================================//
    //======    Here we try to find the events with overlapping jets  =====//
    int n_Overlapping=0; int n_GenJets=0;
    std::vector<int> OverlappingIndex; OverlappingIndex.clear();
    std::vector<int> OverlappingGenIndex; OverlappingGenIndex.clear();

    //Try to count how many (in total) overlapping hadrons I have in each event
    if (testbHadcommonJetIndex.size()>1){ // We count the overlapping hadrons only for the cases we have at least 2 Jets
        for(size_t i=0;i!=testbHadcommonJetIndex.size();i++){
            for(size_t j=i+1;j!=testbHadcommonJetIndex.size();j++){
                if(testbHadcommonJetIndex.at(i)==testbHadcommonJetIndex.at(j)){
                    if (OverlappingIndex.size()<1){
                        OverlappingGenIndex.push_back(testbHadcommonJetIndex.at(j));
                        OverlappingIndex.push_back(j);
                    }  
                    else if (OverlappingIndex.size()>0){
                        for (size_t k=0; k!=OverlappingIndex.size(); k++){ 
                            //if (j==OverlappingIndex.at(k)){
                            if (static_cast<int>(j)==OverlappingIndex.at(k)){
                                isInOverlappingJetsIndex = true;
                                break;
                            }
                            //else if (j!=OverlappingIndex.at(k)) isInOverlappingJetsIndex = false;
                            else if (static_cast<int>(j)!=OverlappingIndex.at(k)) isInOverlappingJetsIndex = false; 
                        }// for k
                        for (size_t g=0; g!=OverlappingGenIndex.size(); g++){
                            if (testbHadcommonJetIndex.at(j)==OverlappingGenIndex.at(g)){ 
                                isInOverlappingGenJetsIndex = true;
                                break;
                           }
                           else if (testbHadcommonJetIndex.at(j)!=OverlappingGenIndex.at(g)) isInOverlappingGenJetsIndex = false;
                        }// for g 

              
                        if (!isInOverlappingJetsIndex) OverlappingIndex.push_back(j); 
                        if (!isInOverlappingGenJetsIndex)  OverlappingGenIndex.push_back(testbHadcommonJetIndex.at(j));
                    }//else if (OverlappingIndex.size()>0) 
            
                    if (OverlappingIndex.size()>0){
                        for(size_t n=0; n!=OverlappingIndex.size();n++){
                            //if(i!=OverlappingIndex.at(n)){
                            if(static_cast<int>(i)!=OverlappingIndex.at(n)){
                                isInOverlappingJetsIndex = false;
                            }//if(i!=OverlappingIndex.at(n))
                            //else if(i==OverlappingIndex.at(n)){
                            else if(static_cast<int>(i)==OverlappingIndex.at(n)){
                                isInOverlappingJetsIndex = true;
                                break;
                            }//else if(i==OverlappingIndex.at(n))   
                        }//for n 
                    }//if (OverlappingIndex.size()>0)
            
            
                    if (!isInOverlappingJetsIndex){ 
                        n_Overlapping=OverlappingIndex.size()+OverlappingGenIndex.size();
                        n_GenJets=OverlappingGenIndex.size();
                    }//if (!isInOverlappingJetsIndex){
                }//if (testbHadcommonJetIndex.at(i)==testbHadcommonJetIndex.at(j))
          
            }//for j
        }//for i

        if (n_Overlapping>0) isOverlapping=true;
        else if (n_Overlapping==0) isOverlapping=false;
        
        name = "NOverlappingHadrons";
        m_histogram[name]->Fill(n_Overlapping);
            
        name = "HasOverlappingHadrons";
        if (isOverlapping)  m_histogram[name]->Fill(1);
        if (!isOverlapping) m_histogram[name]->Fill(0);
      
        name = "GenJetsID_VS_NOverlappingHadrons_2D";
        ((TH2D*)m_histogram[name])->Fill(n_Overlapping,n_GenJets);
          
    }//if (testbHadcommonJetIndex.size()>1) 
    //======================================================================//

    //======    Here we try to find the events with overlapping jets  =====//
    int n_Overlapping_Top=0; int n_TopGenJets=0;
    std::vector<int> OverlappingTopJetsIndex; OverlappingTopJetsIndex.clear();
    std::vector<int> OverlappingTopJetsGenIndex; OverlappingTopJetsGenIndex.clear();

    //Try to count how many (in total) overlapping hadrons I have in each event
    if (MyTopJetsIndex.size()>1){ // We count the overlapping hadrons only for the cases we have at least 2 Top Jets
        for(size_t i=0;i!=MyTopJetsIndex.size();i++){
            for(size_t j=i+1;j!=MyTopJetsIndex.size();j++){
                if(MyTopJetsIndex.at(i)==MyTopJetsIndex.at(j)){
                    if (OverlappingTopJetsIndex.size()<1){
                        OverlappingTopJetsGenIndex.push_back(MyTopJetsIndex.at(j));
                        OverlappingTopJetsIndex.push_back(j);
                    }  
                    else if (OverlappingTopJetsIndex.size()>0){
                        for (size_t k=0; k!=OverlappingTopJetsIndex.size(); k++){ 
                            //if (j==OverlappingTopJetsIndex.at(k)){
                            if (static_cast<int>(j)==OverlappingTopJetsIndex.at(k)){
                                isInOverlappingTopJetsIndex = true;
                                break;
                            }
                            //else if (j!=OverlappingTopJetsIndex.at(k)) isInOverlappingTopJetsIndex = false;
                            else if (static_cast<int>(j)!=OverlappingTopJetsIndex.at(k)) isInOverlappingTopJetsIndex = false; 
                        }// for k
                        for (size_t g=0; g!=OverlappingTopJetsGenIndex.size(); g++){
                            if (MyTopJetsIndex.at(j)==OverlappingTopJetsGenIndex.at(g)){ 
                                isInOverlappingTopGenJetsIndex = true;
                                break;
                            }
                            else if (MyTopJetsIndex.at(j)!=OverlappingTopJetsGenIndex.at(g)) isInOverlappingTopGenJetsIndex = false;
                        }// for g 

              
                        if (!isInOverlappingTopJetsIndex) OverlappingTopJetsIndex.push_back(j); 
                        if (!isInOverlappingTopGenJetsIndex)  OverlappingTopJetsGenIndex.push_back(MyTopJetsIndex.at(j));
                    }//else if (OverlappingTopJetsIndex.size()>0) 
            
                    if (OverlappingTopJetsIndex.size()>0){
                        for(size_t n=0; n!=OverlappingTopJetsIndex.size();n++){
                            //if(i!=OverlappingTopJetsIndex.at(n)){
                            if(static_cast<int>(i)!=OverlappingTopJetsIndex.at(n)){
                                isInOverlappingTopJetsIndex = false;
                            }//if(i!=OverlappingTopJetsIndex.at(n))
                            //else if(i==OverlappingTopJetsIndex.at(n)){
                            else if(static_cast<int>(i)==OverlappingTopJetsIndex.at(n)){
                                isInOverlappingTopJetsIndex = true;
                                break;
                            }//else if(i==OverlappingTopJetsIndex.at(n))   
                        }//for n 
                    }//if (OverlappingTopJetsIndex.size()>0)
             
                    if (!isInOverlappingTopJetsIndex){ 
                        n_Overlapping_Top=OverlappingTopJetsIndex.size()+OverlappingTopJetsGenIndex.size();
                        n_TopGenJets=OverlappingTopJetsGenIndex.size();
                    }//if (MyTopJetsIndex.at(i)==MyTopJetsIndex.at(j))
                }//if (!isInOverlappingTopJetsIndex){
            }//for j
        }//for i

        if (n_Overlapping_Top>0) isOverlapping_Top=true;
        else if (n_Overlapping_Top==0) isOverlapping_Top=false;


        name = "NOverlappingHadrons_Top";
        m_histogram[name]->Fill(n_Overlapping_Top);
        
        name = "HasOverlappingHadrons_Top";
        if (isOverlapping_Top)  m_histogram[name]->Fill(1);
        if (!isOverlapping_Top) m_histogram[name]->Fill(0);
      
        name = "GenJetsID_VS_NOverlappingHadrons_Top_2D";
        ((TH2D*)m_histogram[name])->Fill(n_Overlapping_Top,n_TopGenJets);    
    }//if (MyTopJetsIndex.size()>1)      
    

    //Higgs Jets
    //======    Here we try to find the events with overlapping jets  =====//
    int n_Overlapping_Higgs=0; int n_HiggsGenJets=0;
    std::vector<int> OverlappingHiggsJetsIndex; OverlappingHiggsJetsIndex.clear();
    std::vector<int> OverlappingHiggsJetsGenIndex; OverlappingHiggsJetsGenIndex.clear();

    //Try to count how many (in total) overlapping hadrons I have in each event
    if (MyHiggsJetsIndex.size()>1){ // We count the overlapping hadrons only for the cases we have at least 2 Higgs Jets
        for(size_t i=0;i!=MyHiggsJetsIndex.size();i++){
            for(size_t j=i+1;j!=MyHiggsJetsIndex.size();j++){
                if(MyHiggsJetsIndex.at(i)==MyHiggsJetsIndex.at(j)){ 
                    if (OverlappingHiggsJetsIndex.size()<1){
                        OverlappingHiggsJetsGenIndex.push_back(MyHiggsJetsIndex.at(j));
                        OverlappingHiggsJetsIndex.push_back(j);
                    }  
                    else if (OverlappingHiggsJetsIndex.size()>0){
                        for (size_t k=0; k!=OverlappingHiggsJetsIndex.size(); k++){ 
                            //if (j==OverlappingHiggsJetsIndex.at(k)){
                            if (static_cast<int>(j)==OverlappingHiggsJetsIndex.at(k)){
                                isInOverlappingHiggsJetsIndex = true;
                                break;
                            }
                            //else if (j!=OverlappingHiggsJetsIndex.at(k)) isInOverlappingHiggsJetsIndex = false;
                            else if (static_cast<int>(j)!=OverlappingHiggsJetsIndex.at(k)) isInOverlappingHiggsJetsIndex = false; 
                        }// for k
                        for (size_t g=0; g!=OverlappingHiggsJetsGenIndex.size(); g++){
                            if (MyHiggsJetsIndex.at(j)==OverlappingHiggsJetsGenIndex.at(g)){ 
                                isInOverlappingHiggsGenJetsIndex = true;
                                break;
                            }
                            else if (MyHiggsJetsIndex.at(j)!=OverlappingHiggsJetsGenIndex.at(g)) isInOverlappingHiggsGenJetsIndex = false;
                        }// for g 

                        if (!isInOverlappingHiggsJetsIndex) OverlappingHiggsJetsIndex.push_back(j); 
                        if (!isInOverlappingHiggsGenJetsIndex)  OverlappingHiggsJetsGenIndex.push_back(MyHiggsJetsIndex.at(j));
                    }//else if (OverlappingHiggsJetsIndex.size()>0) 
        
                    if (OverlappingHiggsJetsIndex.size()>0){
                        for(size_t n=0; n!=OverlappingHiggsJetsIndex.size();n++){
                            //if(i!=OverlappingHiggsJetsIndex.at(n)){
                            if(static_cast<int>(i)!=OverlappingHiggsJetsIndex.at(n)){
                                isInOverlappingHiggsJetsIndex = false;
                            }//if(i!=OverlappingHiggsJetsIndex.at(n))
                            //else if(i==OverlappingHiggsJetsIndex.at(n)){
                            else if(static_cast<int>(i)==OverlappingHiggsJetsIndex.at(n)){
                                isInOverlappingHiggsJetsIndex = true;
                                break;
                            }//else if(i==OverlappingHiggsJetsIndex.at(n))   
                        }//for n 
                    }//if (OverlappingHiggsJetsIndex.size()>0)
        
                    if (!isInOverlappingHiggsJetsIndex){ 
                        n_Overlapping_Higgs=OverlappingHiggsJetsIndex.size()+OverlappingHiggsJetsGenIndex.size();
                        n_HiggsGenJets=OverlappingHiggsJetsGenIndex.size();
                    }
                }
            }
        }

    
        if (n_Overlapping_Higgs>0) isOverlapping_Higgs=true;
        else if (n_Overlapping_Higgs==0) isOverlapping_Higgs=false;


        name = "NOverlappingHadrons_Higgs";
        m_histogram[name]->Fill(n_Overlapping_Higgs);
      
        name = "HasOverlappingHadrons_Higgs";
        if (isOverlapping_Higgs)  m_histogram[name]->Fill(1);
        if (!isOverlapping_Higgs) m_histogram[name]->Fill(0);
    
        name = "GenJetsID_VS_NOverlappingHadrons_Higgs_2D";
        ((TH2D*)m_histogram[name])->Fill(n_Overlapping_Higgs,n_HiggsGenJets);
    }
          

    //Top+Higgs Jets
    //======    Here we try to find the events with overlapping jets  =====//
    int n_Overlapping_TopHiggs=0; int n_TopHiggsGenJets=0;
    std::vector<int> OverlappingTopHiggsJetsIndex; OverlappingTopHiggsJetsIndex.clear();
    std::vector<int> OverlappingTopHiggsJetsGenIndex; OverlappingTopHiggsJetsGenIndex.clear();

    //Try to count how many (in total) overlapping hadrons I have in each event
    if (MyTopHiggsJetsIndex.size()>3 && MyTopJetsIndex.size()>1 && MyHiggsJetsIndex.size()>1){ // We count the overlapping hadrons only for the cases we have at least 2 Top + at least 2 Higgs Jets
    // if (MyTopHiggsJetsIndex.size()>1){
        for(size_t i=0;i!=MyTopHiggsJetsIndex.size();i++){
            for(size_t j=i+1;j!=MyTopHiggsJetsIndex.size();j++){
                if(MyTopHiggsJetsIndex.at(i)==MyTopHiggsJetsIndex.at(j)){             
                    if (OverlappingTopHiggsJetsIndex.size()<1){
                        OverlappingTopHiggsJetsGenIndex.push_back(MyTopHiggsJetsIndex.at(j));
                        OverlappingTopHiggsJetsIndex.push_back(j);
                    }  
                    else if (OverlappingTopHiggsJetsIndex.size()>0){
                        for (size_t k=0; k!=OverlappingTopHiggsJetsIndex.size(); k++){ 
                            //if (j==OverlappingTopHiggsJetsIndex.at(k)){
                            if (static_cast<int>(j)==OverlappingTopHiggsJetsIndex.at(k)){
                                isInOverlappingTopHiggsJetsIndex = true;
                                break;
                            }
                            //else if (j!=OverlappingTopHiggsJetsIndex.at(k)) isInOverlappingTopHiggsJetsIndex = false;
                            else if (static_cast<int>(j)!=OverlappingTopHiggsJetsIndex.at(k)) isInOverlappingTopHiggsJetsIndex = false; 
                        }
                        for (size_t g=0; g!=OverlappingTopHiggsJetsGenIndex.size(); g++){
                            if (MyTopHiggsJetsIndex.at(j)==OverlappingTopHiggsJetsGenIndex.at(g)){ 
                               isInOverlappingTopHiggsGenJetsIndex = true;
                               break;
                            }
                            else if (MyTopHiggsJetsIndex.at(j)!=OverlappingTopHiggsJetsGenIndex.at(g)) isInOverlappingTopHiggsGenJetsIndex = false;
                        }
              
                        if (!isInOverlappingTopHiggsJetsIndex) OverlappingTopHiggsJetsIndex.push_back(j); 
                        if (!isInOverlappingTopHiggsGenJetsIndex)  OverlappingTopHiggsJetsGenIndex.push_back(MyTopHiggsJetsIndex.at(j));
                    }
            
                    if (OverlappingTopHiggsJetsIndex.size()>0){
                        for(size_t n=0; n!=OverlappingTopHiggsJetsIndex.size();n++){
                            //if(i!=OverlappingTopHiggsJetsIndex.at(n)){
                            if(static_cast<int>(i)!=OverlappingTopHiggsJetsIndex.at(n)){
                                isInOverlappingTopHiggsJetsIndex = false;
                            }
                            //else if(i==OverlappingTopHiggsJetsIndex.at(n)){
                            else if(static_cast<int>(i)==OverlappingTopHiggsJetsIndex.at(n)){
                                isInOverlappingTopHiggsJetsIndex = true;
                                break;
                            }
                        }
                    }
                        
                    if (!isInOverlappingTopHiggsJetsIndex){ 
                        n_Overlapping_TopHiggs=OverlappingTopHiggsJetsIndex.size()+OverlappingTopHiggsJetsGenIndex.size();
                        n_TopHiggsGenJets=OverlappingTopHiggsJetsGenIndex.size();
                    }
                }
            }
        }
     
        if (n_Overlapping_TopHiggs>0) isOverlapping_TopHiggs=true;
        else if (n_Overlapping_TopHiggs==0) isOverlapping_TopHiggs=false;

        name = "NOverlappingHadrons_TopHiggs";
        m_histogram[name]->Fill(n_Overlapping_TopHiggs);
            
        name = "HasOverlappingHadrons_TopHiggs";
        if (isOverlapping_TopHiggs)  m_histogram[name]->Fill(1);
        if (!isOverlapping_TopHiggs) m_histogram[name]->Fill(0);
      
        name = "GenJetsID_VS_NOverlappingHadrons_TopHiggs_2D";
        ((TH2D*)m_histogram[name])->Fill(n_Overlapping_TopHiggs,n_TopHiggsGenJets);
    }


    //Find the leading Jet  
    //   LeadingJetIndx = testbHadcommonJetIndex.at(0);
    //   LeadingFlavourIndx = 0;

    //   if (MyTopJetsIndex.size()>0) LeadingTopJetIndx = MyTopJetsIndex.at(0);
    //   if (MyHiggsJetsIndex.size()>0) LeadingHiggsJetIndx = MyHiggsJetsIndex.at(0);


    if (!isOverlapping_Top && genObjectIndices.uniqueGenTopMatching()){  // we keep only the events with no overlapping jets!!!
        if (MyTopJetsIndex.size()>0) LeadingTopJetIndx = MyTopJetsIndex.at(0);
      
        //FIXME :: When we have ONLY one Top (Higgs) Jet. In this case the leading and subleading jet is the same jet. 
        //Find the Top Subleading Jet
        if (MyTopJetsIndex.size()>1){
            if(MyTopJetsIndex.at(1) != MyTopJetsIndex.at(0)) SubLeadingTopJetIndx = MyTopJetsIndex.at(1);
            else if (MyTopJetsIndex.size()==2 && MyTopJetsIndex.at(1)==MyTopJetsIndex.at(0)){
                std::cout<< "The two hadrons correspond to the one and only Hadronic Jet." <<std::endl;
                TopHasLeadingAndSubleading = false;
                SubLeadingTopJetIndx = MyTopJetsIndex.at(1);
            }

            else if (MyTopJetsIndex.size()>2 && MyTopJetsIndex.at(1)==MyTopJetsIndex.at(0) && MyTopJetsIndex.at(2)!=MyTopJetsIndex.at(1)) SubLeadingTopJetIndx = MyTopJetsIndex.at(2);
            else if (MyTopJetsIndex.size()==3 && MyTopJetsIndex.at(1)==MyTopJetsIndex.at(0)&& MyTopJetsIndex.at(2)==MyTopJetsIndex.at(1)){
                std::cout<< "The three hadrons correspond to the one and only Hadronic Jet." <<std::endl;
                TopHasLeadingAndSubleading = false;
                SubLeadingTopJetIndx = MyTopJetsIndex.at(1);
            }  
            else SubLeadingTopJetIndx = MyTopJetsIndex.at(3);
        }
    }

    //Find the Higgs Subleading Jet
    if (!isOverlapping_Higgs && genObjectIndices.uniqueGenHiggsMatching()){ // we keep only the events with no overlapping jets!!!
        if (MyHiggsJetsIndex.size()>0) LeadingHiggsJetIndx = MyHiggsJetsIndex.at(0);
        if (MyHiggsJetsIndex.size()>1){
            if(MyHiggsJetsIndex.at(1) != MyHiggsJetsIndex.at(0)) SubLeadingHiggsJetIndx = MyHiggsJetsIndex.at(1);
            else if (MyHiggsJetsIndex.size()==2 && MyHiggsJetsIndex.at(1)==MyHiggsJetsIndex.at(0)){
                std::cout<< "The two hadrons correspond to the one and only Hadronic Jet." <<std::endl;
                HiggsHasLeadingAndSubleading = false;
                SubLeadingHiggsJetIndx = MyHiggsJetsIndex.at(1);
            }

            else if (MyHiggsJetsIndex.size()>2 && MyHiggsJetsIndex.at(1)==MyHiggsJetsIndex.at(0) && MyHiggsJetsIndex.at(2)!=MyHiggsJetsIndex.at(1)) SubLeadingHiggsJetIndx = MyHiggsJetsIndex.at(2);
            else if (MyHiggsJetsIndex.size()==3 && MyHiggsJetsIndex.at(1)==MyHiggsJetsIndex.at(0)&& MyHiggsJetsIndex.at(2)==MyHiggsJetsIndex.at(1)){
                std::cout<< "The three hadrons correspond to the one and only Hadronic Jet." <<std::endl;
                HiggsHasLeadingAndSubleading = false;
                SubLeadingHiggsJetIndx = MyHiggsJetsIndex.at(1);
            }  
            else SubLeadingHiggsJetIndx = MyHiggsJetsIndex.at(3);
        }
    }

     
      
    // Here we find the Minimum values for DeltaEta, DeltaPhi and DeltaR
    //  if(MyTopJetsIndex.size()>1 && MyHiggsJetsIndex.size()>1 && MyTopHiggsJetsIndex.size()>3 && HiggsHasLeadingAndSubleading && !isOverlapping_TopHiggs && genObjectIndices.uniqueGenTopMatching() && genObjectIndices.uniqueGenHiggsMatching()){  
    if(MyTopJetsIndex.size()>1 && MyHiggsJetsIndex.size()>1 && HiggsHasLeadingAndSubleading && !isOverlapping_TopHiggs && genObjectIndices.uniqueGenMatching() ){  
        if (MyTopHiggsJetsIndex.size()<4) std::cout<<"MyTopHiggsJetsIndex.size()= "<<MyTopHiggsJetsIndex.size()<<std::endl;
        minR = 999.;
        minEta= 999.;
        minPhi= 999.;

        for(size_t i_jet=0;i_jet!=MyTopHiggsJetsIndex.size();i_jet++){ 
            for (size_t j_jet=i_jet+1;j_jet!=MyTopHiggsJetsIndex.size();j_jet++){
                deltaR = ROOT::Math::VectorUtil::DeltaR(commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(i_jet)),commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(j_jet)));
                if (deltaR==999.) std::cout<<"ERROR!!"<<std::endl;

                //choose smallest R
                if (deltaR>minR) continue;
                if (deltaR<minR) minR = deltaR;
            }

            for (size_t j_jet=i_jet+1;j_jet!=MyTopHiggsJetsIndex.size();j_jet++){
                deltaEta = (commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(i_jet)).Eta(),commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(j_jet)).Eta());
                if (deltaEta!=999.){
                    deltaEta_Abs = std::abs(deltaEta); 
                }
                else if (deltaEta==999.) deltaEta_Abs = deltaEta;
                if (deltaEta_Abs==999.) std::cout<<"ERROR!!"<<std::endl;

               //choose smallest R
               if (deltaEta_Abs>minEta) continue;
               if (deltaEta_Abs<minEta) minEta = deltaEta_Abs;
            }

            for (size_t j_jet=i_jet+1;j_jet!=MyTopHiggsJetsIndex.size();j_jet++){
                deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(i_jet)),commonGenObjects.allGenJets_->at(MyTopHiggsJetsIndex.at(j_jet)));
                if (deltaPhi!=999.){
                    deltaPhi_Abs = std::abs(deltaPhi);
                }
                else if (deltaPhi==999.) deltaPhi_Abs = deltaPhi;
                if (deltaPhi_Abs==999.) std::cout<<"ERROR!!"<<std::endl;
              
                //choose smallest R
                if (deltaPhi_Abs>minPhi) continue;
                if (deltaPhi_Abs<minPhi) minPhi = deltaPhi_Abs;
            }
        }

        
        name = "MinDeltaR";
        m_histogram[name]->Fill(minR, weight);

        name = "MinDeltaEta";
        m_histogram[name]->Fill(minEta, weight);
        
        name = "MinDeltaPhi";
        m_histogram[name]->Fill(minPhi, weight);

        name = "Topleading_TopHiggsJets_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingTopJetIndx).Pt(), weight);
        
        name = "Topleading_TopHiggsJets_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingTopJetIndx).Eta(), weight);

        name = "Topsubleading_TopHiggsJets_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx).Pt(), weight);

        name = "Topsubleading_TopHiggsJets_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx).Eta(), weight);

        name = "Higgsleading_TopHiggsJets_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx).Pt(), weight);

        name = "Higgsleading_TopHiggsJets_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx).Eta(), weight);

        name = "Higgssubleading_TopHiggsJets_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx).Pt(), weight);

        name = "Higgssubleading_TopHiggsJets_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx).Eta(), weight);


    }
     

      
    // Fill the histograms
    //========================================    Top Jets   =========================================//

    if(MyTopJetsIndex.size()>1 && TopHasLeadingAndSubleading && !isOverlapping_Top && genObjectIndices.uniqueGenTopMatching()){       
        deltaR = ROOT::Math::VectorUtil::DeltaR(commonGenObjects.allGenJets_->at(LeadingTopJetIndx),commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx));
        if (deltaR==999.) std::cout<<"ERROR!!"<<std::endl;
        
        deltaEta = (commonGenObjects.allGenJets_->at(LeadingTopJetIndx).Eta()-commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx).Eta());
        if (deltaEta==999.) std::cout<<"ERROR!!"<<std::endl;
        
        deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(commonGenObjects.allGenJets_->at(LeadingTopJetIndx),commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx));
        if (deltaPhi==999.) std::cout<<"ERROR!!"<<std::endl;


        name = "DeltaR_Topjet";
        m_histogram[name]->Fill(deltaR, weight);
        
        name = "DeltaEta_Topjet";
        m_histogram[name]->Fill(deltaEta, weight);

        name = "DeltaPhi_Topjet";
        m_histogram[name]->Fill(deltaPhi, weight);

        name = "leading_Topjet_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingTopJetIndx).Pt(), weight);

        name = "leading_Topjet_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingTopJetIndx).Eta(), weight);

        name = "subleading_Topjet_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx).Pt(), weight);

        name = "subleading_Topjet_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingTopJetIndx).Eta(), weight);
    }

    //=================================     Higgs Jets    ===========================================//
    if (MyHiggsJetsIndex.size()>1 && HiggsHasLeadingAndSubleading && !isOverlapping_Higgs && genObjectIndices.uniqueGenHiggsMatching()){        
         
        deltaR = ROOT::Math::VectorUtil::DeltaR(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx),commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx));
        if (deltaR==999.) std::cout<<"ERROR!!"<<std::endl;
        
        deltaEta = (commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx).Eta()-commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx).Eta());
        if (deltaEta==999.) std::cout<<"ERROR!!"<<std::endl;
        
        deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx),commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx));
        if (deltaPhi==999.) std::cout<<"ERROR!!"<<std::endl;

        
        name = "DeltaR_Higgsjet";
        m_histogram[name]->Fill(deltaR, weight);

        name = "DeltaEta_Higgsjet";
        m_histogram[name]->Fill(deltaEta, weight);

        name = "DeltaPhi_Higgsjet";
        m_histogram[name]->Fill(deltaPhi, weight);
        
        name = "leading_Higgsjet_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx).Pt(), weight);

        name = "leading_Higgsjet_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(LeadingHiggsJetIndx).Eta(), weight);

        name = "subleading_Higgsjet_Pt";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx).Pt(), weight);

        name = "subleading_Higgsjet_eta";
        m_histogram[name]->Fill(commonGenObjects.allGenJets_->at(SubLeadingHiggsJetIndx).Eta(), weight);
    }


     //================================   Gen - Reco Plots  ===================================================//

    if(TopgenIndex.size()>0 && !isOverlapping_Top && genObjectIndices.uniqueGenTopMatching()){
        for(size_t i_get=0;i_get!=TopgenIndex.size();i_get++){
            if (isInMismatchedTopJetsIndex.size()>0){
                if(isInMismatchedTopJetsIndex.at(i_get)) continue;
                double TopgenP  = commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).M();
                double ToprecoP = recoObjects.jets_->at(ToprecoIndex.at(i_get)).M();
          
                name = "GenP_VS_RecoP_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(ToprecoP,TopgenP,weight);
          
                name = "MinDeltaR_RecoGen_Topjet";
                m_histogram[name]->Fill(minRTopRecoGen.at(i_get),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(), minRTopRecoGen.at(i_get),weight);

                name = "RecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtVSRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);
      
                name = "GenEtaVSRecoEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Eta(),weight);

                name = "GenPhiVSRecoPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Eta(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Phi(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(i_get), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);
            }
        }
    }


    if(HiggsgenIndex.size()>0 && !isOverlapping_Higgs && genObjectIndices.uniqueGenHiggsMatching()){
        for(size_t i_get=0;i_get!=HiggsgenIndex.size();i_get++){
            if (isInMismatchedHiggsJetsIndex.size()>0){
                if(isInMismatchedHiggsJetsIndex.at(i_get)) continue;
                double HiggsgenP  = commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).M();
                double HiggsrecoP = recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).M();          
                
                name = "GenP_VS_RecoP_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(HiggsrecoP,HiggsgenP,weight);
                
                name = "MinDeltaR_RecoGen_Higgsjet";
                m_histogram[name]->Fill(minRHiggsRecoGen.at(i_get),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(), minRHiggsRecoGen.at(i_get),weight);

                name = "RecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtVSRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);
          
                name = "GenEtaVSRecoEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Eta(),weight);

                name = "GenPhiVSRecoPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Eta(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Phi(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(i_get), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);
            }
        }
    }


    if(genIndex.size()>0 && !isOverlapping){
        for(size_t i_get=0;i_get!=genIndex.size();i_get++){
            if (isInMismatchedJetsIndex.size()>0){
                if(isInMismatchedJetsIndex.at(i_get)) continue;
                double genP  = commonGenObjects.allGenJets_->at(genIndex.at(i_get)).M();
                double recoP = recoObjects.jets_->at(recoIndex.at(i_get)).M();          
                
                name = "GenP_VS_RecoP_2D";
                ((TH2D*)m_histogram[name])->Fill(recoP,genP,weight);
                
                name = "MinDeltaR_RecoGen";
                m_histogram[name]->Fill(minRRecoGen.at(i_get),weight);

                name = "MinDeltaR_RecoGen_VS_PtGen_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(), minRRecoGen.at(i_get),weight);

                name = "RecoPtOverGenPt";
                m_histogram[name]->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPt";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "GenPtVSRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);
          
                name = "GenEtaVSRecoEta_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Eta(),weight);

                name = "GenPhiVSRecoPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Phi(),weight);

                name = "GenPtVSGenPtMinusRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenEta_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Eta(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsGenPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Phi(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "RecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(i_get), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);
            }
        }
    }

    //Plots for the WRONG Matching
    if(genIndex.size()>0 && !isOverlapping){
        for(size_t i_get=0;i_get!=genIndex.size();i_get++){
            if (isInMismatchedJetsIndex.size()>0){
                if(!isInMismatchedJetsIndex.at(i_get)) continue;
                double genP  = commonGenObjects.allGenJets_->at(genIndex.at(i_get)).M();
                double recoP = recoObjects.jets_->at(recoIndex.at(i_get)).M();  
                
                name = "Mismatched_GenP_VS_RecoP_2D";
                ((TH2D*)m_histogram[name])->Fill(recoP,genP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen";
                m_histogram[name]->Fill(minRRecoGen.at(i_get),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(), minRRecoGen.at(i_get),weight);

                name = "Mismatched_RecoPtOverGenPt";
                m_histogram[name]->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(recoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Eta(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Phi(), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(i_get), recoObjects.jets_->at(recoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_2D";
                ((TH2D*)m_histogram[name])->Fill(minRRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt()-recoObjects.jets_->at(recoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(genIndex.at(i_get)).Pt(),weight);
            }
        }
    } 
    
    if(TopgenIndex.size()>0 && !isOverlapping_Top && genObjectIndices.uniqueGenTopMatching()){
        for(size_t i_get=0;i_get!=TopgenIndex.size();i_get++){
            if (isInMismatchedTopJetsIndex.size()>0){
                if(!isInMismatchedTopJetsIndex.at(i_get)) continue;
                double TopgenP  = commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).M();
                double ToprecoP = recoObjects.jets_->at(ToprecoIndex.at(i_get)).M();
                
                name = "Mismatched_GenP_VS_RecoP_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(ToprecoP,TopgenP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen_Topjet";
                m_histogram[name]->Fill(minRTopRecoGen.at(i_get),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(), minRTopRecoGen.at(i_get),weight);

                name = "Mismatched_RecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt_Topjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(ToprecoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Eta(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Phi(), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(i_get), recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Topjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRTopRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(ToprecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(TopgenIndex.at(i_get)).Pt(),weight);
            }
        }
    }


    if(HiggsgenIndex.size()>0 && !isOverlapping_Higgs && genObjectIndices.uniqueGenHiggsMatching()){
        for(size_t i_get=0;i_get!=HiggsgenIndex.size();i_get++){
            if (isInMismatchedHiggsJetsIndex.size()>0){
                if(!isInMismatchedHiggsJetsIndex.at(i_get)) continue;
                double HiggsgenP  = commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).M();
                double HiggsrecoP = recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).M();  
                
                name = "Mismatched_GenP_VS_RecoP_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(HiggsrecoP,HiggsgenP,weight);
                
                name = "Mismatched_MinDeltaR_RecoGen_Higgsjet";
                m_histogram[name]->Fill(minRHiggsRecoGen.at(i_get),weight);

                name = "Mismatched_MinDeltaR_RecoGen_VS_PtGen_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(), minRHiggsRecoGen.at(i_get),weight);

                name = "Mismatched_RecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPt_Higgsjet";
                m_histogram[name]->Fill((commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtVSRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);
          
                name = "Mismatched_GenEtaVSRecoEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Eta(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Eta(),weight);

                name = "Mismatched_GenPhiVSRecoPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Phi(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Phi(),weight);

                name = "Mismatched_GenPtVSGenPtMinusRecoPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt(), commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPt_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenEta_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Eta(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsGenPhi_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Phi(), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_RecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(i_get), recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt()/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);

                name = "Mismatched_GenPtMinusRecoPtOverGenPtVsMinDeltaR_Higgsjet_2D";
                ((TH2D*)m_histogram[name])->Fill(minRHiggsRecoGen.at(i_get), (commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt()-recoObjects.jets_->at(HiggsrecoIndex.at(i_get)).Pt())/commonGenObjects.allGenJets_->at(HiggsgenIndex.at(i_get)).Pt(),weight);
            }
        }
    }

    
}








