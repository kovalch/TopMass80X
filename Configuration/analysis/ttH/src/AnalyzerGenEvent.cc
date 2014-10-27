#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerGenEvent.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"





AnalyzerGenEvent::AnalyzerGenEvent(const std::vector<TString>& selectionStepsNoCategories,
                                   const std::vector<TString>& stepsForCategories,
                                   const JetCategories* jetCategories):
AnalyzerBase("genEvent_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up gen event analyzer\n";
    std::cout<<"=== Finishing setting up gen event analyzer\n\n";
}



void AnalyzerGenEvent::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    this->bookEventHistos(step, m_histogram);
    
    this->bookMatchingHistos(step, m_histogram);
    
    this->bookBhadronHistos("allB_", step, m_histogram);
    this->bookBhadronHistos("top_", step, m_histogram);
    this->bookBhadronHistos("higgs_", step, m_histogram);
    this->bookBhadronHistos("topHiggs_", step, m_histogram);
    
    this->bookTopOrHiggsHistos("top_", step, m_histogram);
    this->bookTopOrHiggsHistos("higgs_", step, m_histogram);
    
    this->bookTopAndHiggsHistos(step, m_histogram);
}



void AnalyzerGenEvent::bookEventHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "genJet_all_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"All Gen Jet Multiplicity;# jets_{gen}^{all};# events",80,0,80));
    
    name = "genJet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen Jet Multiplicity;# jets_{gen};# events",80,0,80));
    
    name = "genJet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen Jet Pt;jet p_{T};# jets_{gen}",50,0,300));
    
    name = "genJet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen Jet Eta;jet #eta;# jets_{gen}",50,-2.6,2.6));
    
    name = "genJet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Gen Jet Eta;jet #eta;# jets_{gen}",50,-3.2,3.2));
}



void AnalyzerGenEvent::bookMatchingHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "top_matchedBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "bQuark-genJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "genJet-recoJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Matched");
    
    name = "top_unmatchedGenBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Several hadrons");
    
    name = "top_unmatchedRecoBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    
    name = "higgs_matchedBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "bQuark-genJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "genJet-recoJet fail");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Matched");
    
    name = "higgs_unmatchedGenBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Higgs jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Several hadrons");
    
    name = "higgs_unmatchedRecoBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",3,0,3));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Higgs jets overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "2 jets not matched");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "1 jet not matched");
    
    name = "topHiggs_matchedBjet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Top-Higgs gen overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "Unique gen jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Top-Higgs reco overlap");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "Unique reco jets");
}



void AnalyzerGenEvent::bookBhadronHistos(const TString& topOrHiggsName,
                                         const TString& step,
                                         std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = topOrHiggsName + "bhadrons_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"B hadron Multiplicity;# B hadrons;# events",10,0,10));
    
    name = topOrHiggsName + "bjets_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"b Jet Multiplicity;# b jets;# events",10,0,10));
    
    name = topOrHiggsName + "bhadronsPerBjet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";# B hadrons;# jets",10,0,10));
    
    name = topOrHiggsName + "bhadronsVsBjets_multiplicity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,";# B hadrons;# b jets",10,0,10,10,0,10));
    
    
}



void AnalyzerGenEvent::bookTopOrHiggsHistos(const TString& topOrHiggsName,
                                            const TString& step,
                                            std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = topOrHiggsName + "jet1st_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading jet;p_{T} [GeV];Entries",50,0,300));

    name = topOrHiggsName + "jet1st_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading jet;#eta;Entries",50,-2.6,2.6));

    name = topOrHiggsName + "jet2nd_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading jet;p_{T} [GeV];Entries",50,0,300));
    
    name = topOrHiggsName + "jet2nd_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading jet;#eta;Entries",50,-2.6,2.6));


    name = topOrHiggsName + "dijet_deltaEta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Jets;#Delta#eta_{jj};Entries",50,-5,5));
    
    name = topOrHiggsName + "dijet_deltaPhi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Jets;#Delta#phi_{jj};Entries",40,-4,4));

    name = topOrHiggsName + "dijet_deltaR";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Jets;#DeltaR_{jj};Entries",50,0,5));
}



void AnalyzerGenEvent::bookTopAndHiggsHistos(const TString& step,
                                             std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "topHiggs_dijet_deltaRMin";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs Jets;Min_{#DeltaR_{jj}};Entries",260,0,2));

    name = "topHiggs_dijet_deltaEtaMin";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#eta_{jj}};Entries",40,0,4));
    
    name = "topHiggs_dijet_deltaPhiMin";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#phi_{jj}};Entries",32,0,3));
    
    
    name = "topHiggs_topJet1st_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet in Top+Higgs Jets;p_{T} [GeV];Entries",50,0,300));
    
    name = "topHiggs_topJet1st_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Top jet in Top+Higgs Jets;#eta;Entries",50,-2.6,2.6));

    name = "topHiggs_topJet2nd_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet in Top+Higgs Jets;p_{T} [GeV];Entries",50,0,300));

    name = "topHiggs_topJet2nd_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Top jet in Top+Higgs Jets;#eta;Entries",50,-2.6,2.6));
    
    
    name = "topHiggs_higgsJet1st_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet in Top+Higgs Jets;p_{T} [GeV];Entries",50,0,300));

    name = "topHiggs_higgsJet1st_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Leading Higgs jet in Top+Higgs Jets;#eta;Entries",50,-2.6,2.6));

    name = "topHiggs_higgsJet2nd_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subeading Higgs jet in Top+Higgs Jets;p_{T} [GeV];Entries",50,0,300));
    
    name = "topHiggs_higgsJet2nd_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Subleading Higgs jet in Top+Higgs Jets;#eta;Entries",50,-2.6,2.6));
}



void AnalyzerGenEvent::fillHistos(const EventMetadata&,
                                  const RecoObjects&, const CommonGenObjects&,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                  const KinematicReconstructionSolutions&,
                                  const tth::RecoObjectIndices&, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                  const double& weight, const TString&,
                                  std::map<TString, TH1*>& m_histogram)
{
    if(!topGenObjects.valuesSet_) return;
    
    this->fillEventHistos(topGenObjects, genObjectIndices, weight, m_histogram);
    
    this->fillMatchingHistos(genObjectIndices, weight, m_histogram);
    
    this->fillBhadronHistos("allB_", genObjectIndices, weight, m_histogram);
    this->fillBhadronHistos("top_", genObjectIndices, weight, m_histogram);
    this->fillBhadronHistos("higgs_", genObjectIndices, weight, m_histogram);
    this->fillBhadronHistos("topHiggs_", genObjectIndices, weight, m_histogram);
    
    if(genObjectIndices.uniqueGenTopMatching()){
        this->fillTopOrHiggsHistos("top_", topGenObjects, genObjectIndices, weight, m_histogram);
    }
    
    if(genObjectIndices.uniqueGenHiggsMatching()){
        this->fillTopOrHiggsHistos("higgs_", topGenObjects, genObjectIndices, weight, m_histogram);
    }
    
    if(genObjectIndices.uniqueGenMatching()){
        this->fillTopAndHiggsHistos(topGenObjects, genObjectIndices, weight, m_histogram);
    }
}


void AnalyzerGenEvent::fillEventHistos(const TopGenObjects& topGenObjects,
                                       const tth::GenObjectIndices& genObjectIndices,
                                       const double& weight, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    const VLV& allGenJets = (topGenObjects.valuesSet_) ? *topGenObjects.allGenJets_ : VLV(0);
    const std::vector<int> genJetIndices = genObjectIndices.genJetIndices_;
    
    const int nGenJets_all = allGenJets.size();
    const int nGenJets = genJetIndices.size();
    
    
    name = "genJet_all_multiplicity";
    m_histogram.at(name)->Fill(nGenJets_all, weight);
    
    name = "genJet_multiplicity";
    m_histogram.at(name)->Fill(nGenJets, weight);
    
    for(const int jetId : genJetIndices) {
        name = "genJet_pt";
        m_histogram.at(name)->Fill(allGenJets.at(jetId).Pt());
        
        name = "genJet_eta";
        m_histogram.at(name)->Fill(allGenJets.at(jetId).Eta());
        
        name = "genJet_phi";
        m_histogram.at(name)->Fill(allGenJets.at(jetId).Phi());
    }
    
}



void AnalyzerGenEvent::fillMatchingHistos(const tth::GenObjectIndices& genObjectIndices,
                                          const double& weight, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "top_matchedBjet";
    if(!genObjectIndices.uniqueGenTopMatching()) m_histogram.at(name)->Fill(0., weight);
    else if(!genObjectIndices.uniqueRecoTopMatching()) m_histogram.at(name)->Fill(1., weight);
    else m_histogram.at(name)->Fill(2., weight);

    name = "top_unmatchedGenBjet";
    if(genObjectIndices.genBjetFromTopIndex_==-2 || genObjectIndices.genAntiBjetFromTopIndex_==-2) m_histogram.at(name)->Fill(3., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==-1 && genObjectIndices.genAntiBjetFromTopIndex_==-1) m_histogram.at(name)->Fill(1., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==-1 || genObjectIndices.genAntiBjetFromTopIndex_==-1) m_histogram.at(name)->Fill(2., weight);
    else if(genObjectIndices.genBjetFromTopIndex_==genObjectIndices.genAntiBjetFromTopIndex_) m_histogram.at(name)->Fill(0., weight);

    name = "top_unmatchedRecoBjet";
    if(genObjectIndices.uniqueGenTopMatching()){
        if(genObjectIndices.recoBjetFromTopIndex_<0 && genObjectIndices.recoAntiBjetFromTopIndex_<0) m_histogram.at(name)->Fill(1., weight);
        else if(genObjectIndices.recoBjetFromTopIndex_<0 || genObjectIndices.recoAntiBjetFromTopIndex_<0) m_histogram.at(name)->Fill(2., weight);
        else if(genObjectIndices.recoBjetFromTopIndex_==genObjectIndices.recoAntiBjetFromTopIndex_) m_histogram.at(name)->Fill(0., weight);
    }

    name = "higgs_matchedBjet";
    if(!genObjectIndices.uniqueGenHiggsMatching()) m_histogram.at(name)->Fill(0., weight);
    else if(!genObjectIndices.uniqueRecoHiggsMatching()) m_histogram.at(name)->Fill(1., weight);
    else m_histogram.at(name)->Fill(2., weight);

    name = "higgs_unmatchedGenBjet";
    if(genObjectIndices.genBjetFromHiggsIndex_==-2 || genObjectIndices.genAntiBjetFromHiggsIndex_==-2) m_histogram.at(name)->Fill(3., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==-1 && genObjectIndices.genAntiBjetFromHiggsIndex_==-1) m_histogram.at(name)->Fill(1., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==-1 || genObjectIndices.genAntiBjetFromHiggsIndex_==-1) m_histogram.at(name)->Fill(2., weight);
    else if(genObjectIndices.genBjetFromHiggsIndex_==genObjectIndices.genAntiBjetFromHiggsIndex_) m_histogram.at(name)->Fill(0., weight);

    name = "higgs_unmatchedRecoBjet";
    if(genObjectIndices.uniqueGenHiggsMatching()){
        if(genObjectIndices.recoBjetFromHiggsIndex_<0 && genObjectIndices.recoAntiBjetFromHiggsIndex_<0) m_histogram.at(name)->Fill(1., weight);
        else if(genObjectIndices.recoBjetFromHiggsIndex_<0 || genObjectIndices.recoAntiBjetFromHiggsIndex_<0) m_histogram.at(name)->Fill(2., weight);
        else if(genObjectIndices.recoBjetFromHiggsIndex_==genObjectIndices.recoAntiBjetFromHiggsIndex_) m_histogram.at(name)->Fill(0., weight);
    }

    name = "topHiggs_matchedBjet";
    if(genObjectIndices.uniqueGenMatching()){
        m_histogram.at(name)->Fill(1., weight);
        if(genObjectIndices.uniqueRecoMatching()) m_histogram.at(name)->Fill(3., weight);
        else if(genObjectIndices.uniqueRecoTopMatching() && genObjectIndices.uniqueRecoHiggsMatching()) m_histogram.at(name)->Fill(2., weight);
    }
    else if(genObjectIndices.uniqueGenTopMatching() && genObjectIndices.uniqueGenHiggsMatching()) m_histogram.at(name)->Fill(0., weight);
}



void AnalyzerGenEvent::fillBhadronHistos(const TString& topOrHiggsName,
                                         const tth::GenObjectIndices& genObjectIndices,
                                         const double& weight,
                                         std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Access all b jets for requested selection, avoiding double-counting
    std::vector<int> bjetIndices;
    if(topOrHiggsName == "allB_"){
        bjetIndices = common::initialiseIndices(genObjectIndices.allGenBjetIndices_);
    }
    else if(topOrHiggsName == "top_"){
        if(genObjectIndices.genBjetFromTopIndex_ >= 0) bjetIndices.push_back(genObjectIndices.genBjetFromTopIndex_);
        if(genObjectIndices.genAntiBjetFromTopIndex_>=0 &&
           std::find(bjetIndices.begin(), bjetIndices.end(), genObjectIndices.genAntiBjetFromTopIndex_)!=bjetIndices.end())
            bjetIndices.push_back(genObjectIndices.genAntiBjetFromTopIndex_);
    }
    else if(topOrHiggsName == "higgs_"){
        if(genObjectIndices.genBjetFromHiggsIndex_ >= 0) bjetIndices.push_back(genObjectIndices.genBjetFromHiggsIndex_);
        if(genObjectIndices.genAntiBjetFromHiggsIndex_>=0 &&
           std::find(bjetIndices.begin(), bjetIndices.end(), genObjectIndices.genAntiBjetFromTopIndex_)!=bjetIndices.end())
            bjetIndices.push_back(genObjectIndices.genAntiBjetFromHiggsIndex_);
    }
    else if(topOrHiggsName == "topHiggs_"){
        if(genObjectIndices.genBjetFromTopIndex_ >= 0) bjetIndices.push_back(genObjectIndices.genBjetFromTopIndex_);
        if(genObjectIndices.genAntiBjetFromTopIndex_>=0 &&
           std::find(bjetIndices.begin(), bjetIndices.end(), genObjectIndices.genAntiBjetFromTopIndex_)!=bjetIndices.end())
            bjetIndices.push_back(genObjectIndices.genAntiBjetFromTopIndex_);
        if(genObjectIndices.genBjetFromHiggsIndex_>=0 &&
           std::find(bjetIndices.begin(), bjetIndices.end(), genObjectIndices.genAntiBjetFromTopIndex_)!=bjetIndices.end())
            bjetIndices.push_back(genObjectIndices.genBjetFromHiggsIndex_);
        if(genObjectIndices.genAntiBjetFromHiggsIndex_>=0 &&
           std::find(bjetIndices.begin(), bjetIndices.end(), genObjectIndices.genAntiBjetFromTopIndex_)!=bjetIndices.end())
            bjetIndices.push_back(genObjectIndices.genAntiBjetFromHiggsIndex_);
    }
    else{
        std::cerr<<"ERROR in AnalyzerGenEvent::fillBhadronHistos()! Identifier topOrHiggsName is not allowed: "
                 <<topOrHiggsName<<"\n...break\n"<<std::endl;
        exit(87);
    }
    const int nBjets = bjetIndices.size();
    
    // Fill histograms
    name = topOrHiggsName + "bjets_multiplicity";
    m_histogram.at(name)->Fill(nBjets, weight);
    
    int nBhadrons = 0;
    for(const int bjetIndex : bjetIndices){
        const auto& bhadronIndices = genObjectIndices.genJetBhadronIndices_.at(bjetIndex);
        nBhadrons += bhadronIndices.size();
        
        name = topOrHiggsName + "bhadronsPerBjet_multiplicity";
        m_histogram.at(name)->Fill(bhadronIndices.size(), weight);
    }
    
    name = topOrHiggsName + "bhadrons_multiplicity";
    m_histogram.at(name)->Fill(nBhadrons, weight);
    
    name = topOrHiggsName + "bhadronsVsBjets_multiplicity";
    ((TH2*)m_histogram.at(name))->Fill(nBhadrons, nBjets, weight);
}



void AnalyzerGenEvent::fillTopOrHiggsHistos(const TString& topOrHiggsName,
                                            const TopGenObjects& topGenObjects,
                                            const tth::GenObjectIndices& genObjectIndices,
                                            const double& weight,
                                            std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Access the two requested jets, ordered by pt
    int leadingJetIndex(-1);
    int subleadingJetIndex(-1);
    if(topOrHiggsName == "top_"){
        leadingJetIndex = genObjectIndices.genBjetFromTopIndex_;
        subleadingJetIndex = genObjectIndices.genAntiBjetFromTopIndex_;
    }
    else if(topOrHiggsName == "higgs_"){
        leadingJetIndex = genObjectIndices.genBjetFromHiggsIndex_;
        subleadingJetIndex = genObjectIndices.genAntiBjetFromHiggsIndex_;
    }
    else{
        std::cout<<"ERROR in AnalyzerGenEvent::fillTopOrHiggsHistos()! The given topOrHiggsName is not allowed: "
                 <<topOrHiggsName<<"\n...break\n"<<std::endl;
        exit(98);
    }
    common::orderIndices(leadingJetIndex, subleadingJetIndex, *topGenObjects.allGenJets_, common::LVpt);
    const LV& leadingJet(topGenObjects.allGenJets_->at(leadingJetIndex));
    const LV& subleadingJet(topGenObjects.allGenJets_->at(subleadingJetIndex));
    
    // Calculate quantities
    const double deltaEta = leadingJet.eta() - subleadingJet.eta();
    const double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(leadingJet, subleadingJet);
    const double deltaR = ROOT::Math::VectorUtil::DeltaR(leadingJet, subleadingJet);
    
    // Fill histograms
    name = topOrHiggsName + "jet1st_pt";
    m_histogram.at(name)->Fill(leadingJet.pt(), weight);

    name = topOrHiggsName + "jet1st_eta";
    m_histogram.at(name)->Fill(leadingJet.eta(), weight);

    name = topOrHiggsName + "jet2nd_pt";
    m_histogram.at(name)->Fill(subleadingJet.pt(), weight);
    
    name = topOrHiggsName + "jet2nd_eta";
    m_histogram.at(name)->Fill(subleadingJet.eta(), weight);
    
    name = topOrHiggsName + "dijet_deltaEta";
    m_histogram.at(name)->Fill(deltaEta, weight);
    
    name = topOrHiggsName + "dijet_deltaPhi";
    m_histogram.at(name)->Fill(deltaPhi, weight);

    name = topOrHiggsName + "dijet_deltaR";
    m_histogram.at(name)->Fill(deltaR, weight);
}



void AnalyzerGenEvent::fillTopAndHiggsHistos(const TopGenObjects& topGenObjects,
                                             const tth::GenObjectIndices& genObjectIndices,
                                             const double& weight,
                                             std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Access the top jets, ordered by pt
    int topLeadingJetIndex(genObjectIndices.genBjetFromTopIndex_);
    int topSubleadingJetIndex(genObjectIndices.genAntiBjetFromTopIndex_);
    common::orderIndices(topLeadingJetIndex, topSubleadingJetIndex, *topGenObjects.allGenJets_, common::LVpt);
    const LV& topLeadingJet(topGenObjects.allGenJets_->at(topLeadingJetIndex));
    const LV& topSubleadingJet(topGenObjects.allGenJets_->at(topSubleadingJetIndex));
    
    // Access the Higgs jets, ordered by pt
    int higgsLeadingJetIndex(genObjectIndices.genBjetFromHiggsIndex_);
    int higgsSubleadingJetIndex(genObjectIndices.genAntiBjetFromHiggsIndex_);
    common::orderIndices(higgsLeadingJetIndex, higgsSubleadingJetIndex, *topGenObjects.allGenJets_, common::LVpt);
    const LV& higgsLeadingJet(topGenObjects.allGenJets_->at(higgsLeadingJetIndex));
    const LV& higgsSubleadingJet(topGenObjects.allGenJets_->at(higgsSubleadingJetIndex));
    
    // Calculate quantities
    double minDeltaEta = 999.;
    double minDeltaPhi = 999.;
    double minDeltaR = 999.;
    const VLV& v_jet = {topLeadingJet, topSubleadingJet, higgsLeadingJet, higgsSubleadingJet};
    for(VLV::const_iterator i_jet = v_jet.begin(); i_jet != v_jet.end(); ++i_jet){
        for(VLV::const_iterator j_jet = i_jet+1; j_jet != v_jet.end(); ++j_jet){
            const double deltaEtaAbs = std::abs(i_jet->eta() - j_jet->eta());
            const double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(*i_jet, *j_jet);
            const double deltaR = ROOT::Math::VectorUtil::DeltaR(*i_jet, *j_jet);
            if(deltaEtaAbs < minDeltaEta) minDeltaEta = deltaEtaAbs;
            if(deltaPhi < minDeltaPhi) minDeltaPhi = deltaPhi;
            if(deltaR < minDeltaR) minDeltaR = deltaR;
        }
    }
    
    // Fill histograms
    name = "topHiggs_dijet_deltaEtaMin";
    m_histogram.at(name)->Fill(minDeltaEta, weight);
    
    name = "topHiggs_dijet_deltaPhiMin";
    m_histogram.at(name)->Fill(minDeltaPhi, weight);
    
    name = "topHiggs_dijet_deltaRMin";
    m_histogram.at(name)->Fill(minDeltaR, weight);

    
    name = "topHiggs_topJet1st_pt";
    m_histogram.at(name)->Fill(topLeadingJet.pt(), weight);
    
    name = "topHiggs_topJet1st_eta";
    m_histogram.at(name)->Fill(topLeadingJet.eta(), weight);

    name = "topHiggs_topJet2nd_pt";
    m_histogram.at(name)->Fill(topSubleadingJet.pt(), weight);

    name = "topHiggs_topJet2nd_eta";
    m_histogram.at(name)->Fill(topSubleadingJet.eta(), weight);
    
    
    name = "topHiggs_higgsJet1st_pt";
    m_histogram.at(name)->Fill(higgsLeadingJet.pt(), weight);

    name = "topHiggs_higgsJet1st_eta";
    m_histogram.at(name)->Fill(higgsLeadingJet.eta(), weight);

    name = "topHiggs_higgsJet2nd_pt";
    m_histogram.at(name)->Fill(higgsSubleadingJet.pt(), weight);
    
    name = "topHiggs_higgsJet2nd_eta";
    m_histogram.at(name)->Fill(higgsSubleadingJet.eta(), weight);
}







