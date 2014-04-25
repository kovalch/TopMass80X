#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerKinReco.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerKinReco::AnalyzerKinReco(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("KinReco_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerKinReco::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
        name = "nRecoEvt_vs_JetMult";
        m_histogram[name]    = store(new TH1D ( prefix_+name+step, "N reco events vs Jet Multiplicity", 10, -0.5, 9.5 ));
        
        name = "nRecoEvt_vs_LepEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "lead-LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_LepEta2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "nLead-LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_LeppT";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "lead-LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_LeppT2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "nLead-LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        
        name = "nRecoEvt_vs_JetEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_JetEta2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_JetpT";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_JetpT2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        
        name = "nRecoEvt_vs_LeptonEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_AntiLeptonEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "AntiLeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_LeptonpT";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_AntiLeptonpT";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "AntiLeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        
        name = "nRecoEvt_vs_MET";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "MET - Kin. Reco. behaviour;MET;eff.", 10, 0, 400 ));
        
        name = "nRecoEvt_Eff";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "", 1, 0, 2 ));
        

        name = "signalTopEvents_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N top events vs Jet Multiplicity", 10, -0.5, 9.5 ));
        name = "MatchedJets_vs_JetMult";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "N events with 2 jets from t#bar{t}", 10, -0.5, 9.5 ));
    
        name = "nSolTtJets_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N sol with 2 jets from t#bar{t}", 10, -0.5, 9.5 ));
        name = "nSolCorrSignJets_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N sol with correct sign jets", 10, -0.5, 9.5 ));
}



void AnalyzerKinReco::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                      const TopGenObjects& topGenObjects,
                                      const KinRecoObjects& kinRecoObjects,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
    
    double weight7 = (recoLevelWeights.weightNoPileup_)*(genLevelWeights.weightPileup_);
        weight7 *= recoLevelWeights.weightBtagSF_;
    
    int Njets = recoObjectIndices.jetIndices_.size();
    
        m_histogram["nRecoEvt_vs_JetMult"]->Fill(Njets,weight7);
        
        m_histogram["nRecoEvt_vs_LepEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leadingLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LepEta2"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.nLeadingLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LeppT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leadingLeptonIndex_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_LeppT2"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.nLeadingLeptonIndex_).Pt(),weight7);
        
        m_histogram["nRecoEvt_vs_LeptonEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_AntiLeptonEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.antiLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LeptonpT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leptonIndex_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_AntiLeptonpT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.antiLeptonIndex_).Pt(),weight7);
        
        m_histogram["nRecoEvt_vs_MET"]->Fill( (*recoObjects.met_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_JetEta"]->Fill((*recoObjects.jets_).at(0).Eta(),weight7);
        m_histogram["nRecoEvt_vs_JetEta2"]->Fill((*recoObjects.jets_).at(1).Eta(),weight7);
        m_histogram["nRecoEvt_vs_JetpT"]->Fill((*recoObjects.jets_).at(0).Pt(),weight7);
        m_histogram["nRecoEvt_vs_JetpT2"]->Fill((*recoObjects.jets_).at(1).Pt(),weight7);
        
        m_histogram["nRecoEvt_Eff"]->Fill(1,weight7);
        
        
        
    if(topGenObjects.valuesSet_ && kinRecoObjects.valuesSet_){
        
        m_histogram["signalTopEvents_vs_JetMult"]->Fill(Njets,weight);
        
     
        int bjet_index=-1;
        int bbarjet_index=-1;
        int b_matched_jetIndex=-1,bbar_matched_jetIndex=-1;
        
            for(int i=0;i<(int)((*topGenObjects.genBHadFlavour_).size());i++)
            {
                if(((*topGenObjects.genBHadFlavour_).at(i))==6)bjet_index =  (*topGenObjects.genBHadJetIndex_).at(((*topGenObjects.genBHadIndex_).at(i)));
                if(((*topGenObjects.genBHadFlavour_).at(i))==-6)bbarjet_index =  (*topGenObjects.genBHadJetIndex_).at(((*topGenObjects.genBHadIndex_).at(i)));
            }
        if(bjet_index>=0&&bbarjet_index>=0&&bjet_index!=bbarjet_index){
            
            TLorentzVector truejetB = common::LVtoTLV((*commonGenObjects.allGenJets_).at(bjet_index));
            TLorentzVector truejetBbar = common::LVtoTLV((*commonGenObjects.allGenJets_).at(bbarjet_index));
            TLorentzVector bjetHyp=common::LVtoTLV((*kinRecoObjects.HypBJet_).at(0));
            TLorentzVector bBarjetHyp=common::LVtoTLV((*kinRecoObjects.HypAntiBJet_).at(0));
            
            double dR_b_b=truejetB.DeltaR(bjetHyp);
            double dR_bar_bar=truejetBbar.DeltaR(bBarjetHyp);
            double dR_b_bar=truejetB.DeltaR(bBarjetHyp);
            double dR_bar_b=truejetBbar.DeltaR(bjetHyp);
            
                if((dR_b_b<0.3&&dR_bar_bar<0.3)||(dR_b_bar<0.3&&dR_bar_b<0.3))
                {
                        m_histogram["nSolTtJets_vs_JetMult"]->Fill(Njets,weight);
                        if((dR_b_b<0.3&&dR_bar_bar<0.3))
                        {
                            m_histogram["nSolCorrSignJets_vs_JetMult"]->Fill(Njets,weight);
                        }
                }
            
         for(const int index : (recoObjectIndices.jetIndices_))
         { 
            double dr=100;
            TLorentzVector recojet = common::LVtoTLV((*recoObjects.jets_).at(index));
            dr=recojet.DeltaR(common::LVtoTLV((*commonGenObjects.allGenJets_).at(bjet_index)));
            if(dr<0.3)b_matched_jetIndex=index;
            dr=recojet.DeltaR(common::LVtoTLV((*commonGenObjects.allGenJets_).at(bbarjet_index)));
            if(dr<0.3)bbar_matched_jetIndex=index;
         }
         if(b_matched_jetIndex>=0 && bbar_matched_jetIndex>=0 && b_matched_jetIndex!=bbar_matched_jetIndex)
         {
            m_histogram["MatchedJets_vs_JetMult"]->Fill(Njets,weight);
         }
        
        }
    }

}








