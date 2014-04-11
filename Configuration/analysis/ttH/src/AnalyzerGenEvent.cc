#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerGenEvent.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








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
    TString name;
       
    int bins_pt           = 50;  int max_pt           = 300;
    int bins_eta          = 50;  float min_eta        = -2.6;  float max_eta          = 2.6;
    int bins_phi          = 40;  int min_phi          = -4;    int max_phi            = 4;
    int bins_DEta         = 50;  int min_DEta         = -5;    int max_DEta           = 5;
    int bins_MinDEta         = 40;  int min_MinDEta   = 0;   int max_MinDEta          = 4;
    int bins_MinDPhi      = 32;  int min_MinDPhi      = 0;  int max_MinDPhi           = 3;
    int bins_DR           = 50;  int min_DR           = 0;    int max_DR              = 5;
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
  
    ///==================================               Top + Higgs jets              =========================================//
    name = "MinDeltaR";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs Jets;Min_{#DeltaR_{jj}};Entries",bins_MinDR,min_MinDR,max_MinDR));

    name = "MinDeltaEta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#eta_{jj}};Entries",bins_MinDEta,min_MinDEta,max_MinDEta));
    
    name = "MinDeltaPhi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top+Higgs jets;Min_{#Delta#phi_{jj}};Entries",bins_MinDPhi,min_MinDPhi,max_MinDPhi));
    

    ///==================================               Top jets              =========================================//
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


   ///==================================               Higgs Jets              =========================================//
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
        
}



void AnalyzerGenEvent::fillHistos(const RecoObjects&, const CommonGenObjects& commonGenObjects,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                  const KinRecoObjects&,
                                  const tth::RecoObjectIndices&, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                  const double& weight, const TString&,
                                  std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
   
    /// Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& genJets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
    std::vector<int> genJetsIndices = common::initialiseIndices(genJets);
    common::orderIndices(genJetsIndices, genJets, common::LVpt);

    std::vector<int> bHadGenJetIndex; 

    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);

    double deltaR=999.; double deltaEta=999.; double deltaPhi=999.;
                        double deltaEta_Abs=999.; double deltaPhi_Abs=999.9; 
    double minR = 999.; double minEta = 999.; double minPhi = 999.;
        
    ///Loop over the GenJets and choose the jets that are also included in the Hadronic List (bHadJetIndex) 
    if (genJetsIndices.size()>0){
        for(const int Jetindex : genJetsIndices){ 
            for (const int Hadindex : bHadJetIndex){
                if (Jetindex!=Hadindex) continue;
                bHadGenJetIndex.push_back(Hadindex); 
            }
        }
    }
        
    ///============================     Find the Top and Higgs Jets  (Gen Level)  ================================================//
    std::vector<int> topGenJetsIndex;       topGenJetsIndex.clear();
    std::vector<int> higgsGenJetsIndex;     higgsGenJetsIndex.clear();
    std::vector<int> topHiggsGenJetsIndex;  topHiggsGenJetsIndex.clear();


    //Find the Top Jets
    if (genObjectIndices.genBjetFromTopIndex_>=0)     topGenJetsIndex.push_back(genObjectIndices.genBjetFromTopIndex_);   
    if (genObjectIndices.genAntiBjetFromTopIndex_>=0) topGenJetsIndex.push_back(genObjectIndices.genAntiBjetFromTopIndex_);
    
    //Find the Higgs Jets
    if (genObjectIndices.genBjetFromHiggsIndex_>=0)     higgsGenJetsIndex.push_back(genObjectIndices.genBjetFromHiggsIndex_);    
    if (genObjectIndices.genAntiBjetFromHiggsIndex_>=0) higgsGenJetsIndex.push_back(genObjectIndices.genAntiBjetFromHiggsIndex_);

    //Find Jets with Top and Higgs Jets
    for (size_t i=0;i!=topGenJetsIndex.size();i++){
        topHiggsGenJetsIndex.push_back(topGenJetsIndex.at(i));
    }
    for (size_t i=0;i!=higgsGenJetsIndex.size();i++){
        topHiggsGenJetsIndex.push_back(higgsGenJetsIndex.at(i));
    }

    common::orderIndices(topGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
    common::orderIndices(higgsGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
    common::orderIndices(topHiggsGenJetsIndex, genJets, common::LVpt);  //pt ordered jets 
 

    ///============================       Find the events with overlapping hadrons      =========================================//
    
    if(bHadGenJetIndex.size()>1){ //Count the overlapping hadrons only for the cases with at least 2 Jets
        int nOverlapping = overlappingHadrons( bHadGenJetIndex, true);
        int nGenJets = overlappingHadrons( bHadGenJetIndex, false);
        
        name = "NOverlappingHadrons";
        m_histogram[name]->Fill(nOverlapping);
    
        name = "HasOverlappingHadrons";
        if (nOverlapping>0)  m_histogram[name]->Fill(1);
        if (nOverlapping==0) m_histogram[name]->Fill(0);
    
        name = "GenJetsID_VS_NOverlappingHadrons_2D";
        ((TH2D*)m_histogram[name])->Fill(nOverlapping,nGenJets);
    }
    
    //Top Jets
    if(topGenJetsIndex.size()>1){ //Count the overlapping hadrons only for the cases with at least 2 Top Jets
        int nOverlappingTop = overlappingHadrons( topGenJetsIndex, true);
        int nTopGenJets = overlappingHadrons( topGenJetsIndex, false);
        
        name = "NOverlappingHadrons_Top";
        m_histogram[name]->Fill(nOverlappingTop);
    
        name = "HasOverlappingHadrons_Top";
        if (nOverlappingTop>0)  m_histogram[name]->Fill(1);
        if (nOverlappingTop==0) m_histogram[name]->Fill(0);
    
        name = "GenJetsID_VS_NOverlappingHadrons_Top_2D";
        ((TH2D*)m_histogram[name])->Fill(nOverlappingTop,nTopGenJets);
    }
    //Higgs Jets
    if(higgsGenJetsIndex.size()>1){ //Count the overlapping hadrons only for the cases with at least 2 Higgs Jets
        int nOverlappingHiggs = overlappingHadrons( higgsGenJetsIndex, true);
        int nHiggsGenJets = overlappingHadrons( higgsGenJetsIndex, false);
        name = "NOverlappingHadrons_Higgs";
        m_histogram[name]->Fill(nOverlappingHiggs);
    
        name = "HasOverlappingHadrons_Higgs";
        if (nOverlappingHiggs>0)  m_histogram[name]->Fill(1);
        if (nOverlappingHiggs==0) m_histogram[name]->Fill(0);
    
        name = "GenJetsID_VS_NOverlappingHadrons_Higgs_2D";
        ((TH2D*)m_histogram[name])->Fill(nOverlappingHiggs,nHiggsGenJets);
    }
    
    //Top+Higgs Jets
    if(topHiggsGenJetsIndex.size()>3 && topGenJetsIndex.size()>1 && higgsGenJetsIndex.size()>1){ //Count the overlapping hadrons only for the cases with at least 2 Top + at least 2 Higgs Jets
        int nOverlappingTopHiggs = overlappingHadrons( topHiggsGenJetsIndex, true);
        int nTopHiggsGenJets = overlappingHadrons( topHiggsGenJetsIndex, false);
        name = "NOverlappingHadrons_TopHiggs";
        m_histogram[name]->Fill(nOverlappingTopHiggs);
    
        name = "HasOverlappingHadrons_TopHiggs";
        if (nOverlappingTopHiggs>0)  m_histogram[name]->Fill(1);
        if (nOverlappingTopHiggs==0) m_histogram[name]->Fill(0);
    
        name = "GenJetsID_VS_NOverlappingHadrons_TopHiggs_2D";
        ((TH2D*)m_histogram[name])->Fill(nOverlappingTopHiggs,nTopHiggsGenJets);
    }    
    
    ///===================================               Find the leading Jet            =================================//  
    //Top Jets
    int LeadingTopJetIndx =-1;
    int SubLeadingTopJetIndx =-1;
    
    //Higgs Jets
    int LeadingHiggsJetIndx =-1;
    int SubLeadingHiggsJetIndx =-1;
        
    if (genObjectIndices.uniqueGenTopMatching()){
        LeadingTopJetIndx = topGenJetsIndex.at(0);
        SubLeadingTopJetIndx = topGenJetsIndex.at(1);
    }
    
    if (genObjectIndices.uniqueGenHiggsMatching()){
        LeadingHiggsJetIndx = higgsGenJetsIndex.at(0);
        SubLeadingHiggsJetIndx = higgsGenJetsIndex.at(1);
    }

 
 ///============================  Find the Minimum values for DeltaEta, DeltaPhi and DeltaR ===============================//
    if(genObjectIndices.uniqueGenMatching() ){  
        if (topHiggsGenJetsIndex.size()<4) std::cout<<"topHiggsGenJetsIndex.size()= "<<topHiggsGenJetsIndex.size()<<std::endl;
        minR = 999.;
        minEta= 999.;
        minPhi= 999.;

        for(size_t i_jet=0;i_jet!=topHiggsGenJetsIndex.size();i_jet++){ 
            for (size_t j_jet=i_jet+1;j_jet!=topHiggsGenJetsIndex.size();j_jet++){
                deltaR = ROOT::Math::VectorUtil::DeltaR(commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(i_jet)),commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(j_jet)));
                if (deltaR==999.) std::cout<<"ERROR!!"<<std::endl;

                //choose smallest R
                if (deltaR>minR) continue;
                if (deltaR<minR) minR = deltaR;
            }

            for (size_t j_jet=i_jet+1;j_jet!=topHiggsGenJetsIndex.size();j_jet++){
                deltaEta = (commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(i_jet)).Eta(),commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(j_jet)).Eta());
                if (deltaEta!=999.){
                    deltaEta_Abs = std::abs(deltaEta); 
                }
                else if (deltaEta==999.) deltaEta_Abs = deltaEta;
                if (deltaEta_Abs==999.) std::cout<<"ERROR!!"<<std::endl;

               //choose smallest R
               if (deltaEta_Abs>minEta) continue;
               if (deltaEta_Abs<minEta) minEta = deltaEta_Abs;
            }

            for (size_t j_jet=i_jet+1;j_jet!=topHiggsGenJetsIndex.size();j_jet++){
                deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(i_jet)),commonGenObjects.allGenJets_->at(topHiggsGenJetsIndex.at(j_jet)));
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
     

    ///========================================            Top Jets          =========================================//
    if(genObjectIndices.uniqueGenTopMatching()){       
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

    ///=======================================           Higgs Jets        ===========================================//
    if (genObjectIndices.uniqueGenHiggsMatching()){         
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
   
}


int AnalyzerGenEvent::overlappingHadrons( const std::vector<int>& genIndex, const bool& returnNOverlapping)
{

    bool isInOverlappingJetsIndex = false; bool isInOverlappingGenJetsIndex = false;
    int n_Overlapping=0; int n_GenJets=0; 
    std::vector<int> OverlappingIndex; OverlappingIndex.clear();
    std::vector<int> OverlappingGenIndex; OverlappingGenIndex.clear();

    int returnValue = 0;
    
    ///Count how many (in total) overlapping hadrons there are in each event
    for(size_t i=0;i!=genIndex.size();i++){
        for(size_t j=i+1;j!=genIndex.size();j++){
            if(genIndex.at(i)==genIndex.at(j)){
                if (OverlappingIndex.size()<1){
                    OverlappingGenIndex.push_back(genIndex.at(j));
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
                        if (genIndex.at(j)==OverlappingGenIndex.at(g)){ 
                            isInOverlappingGenJetsIndex = true;
                            break;
                        }
                        else if (genIndex.at(j)!=OverlappingGenIndex.at(g)) isInOverlappingGenJetsIndex = false;
                    }
                    
                    
                    if (!isInOverlappingJetsIndex) OverlappingIndex.push_back(j); 
                    if (!isInOverlappingGenJetsIndex)  OverlappingGenIndex.push_back(genIndex.at(j));
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
                    n_GenJets=OverlappingGenIndex.size();
                }
            }
            
        }
    }
    //choose which value to return
    if(returnNOverlapping) returnValue = n_Overlapping;
    else if (!returnNOverlapping) returnValue = n_GenJets;
    return returnValue;

}



