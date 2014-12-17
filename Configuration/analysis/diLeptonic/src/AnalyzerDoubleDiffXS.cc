#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerDoubleDiffXS.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "TLorentzVector.h"

void setBinning(const TString& name,int& nbinX,int& nbinY, std::vector<Double_t>& x_binsArr, std::vector<Double_t>& y_binsArr,int& ngBins)
{
        TString lineTString;
        
    ifstream myfile("HistoList_2d");
         std::string line;
         
         if (myfile.is_open())
     {
        while ( myfile.good() )
        {
            std::getline (myfile,line);
             //skip empty lines and/or comments
             if (line.size() == 0 || line[0] == '#') continue;
            lineTString=line;
            if(lineTString.Contains(name))
            {
                std::stringstream linestream(line);
                std::string trash;
        
                linestream >> trash; //name
                linestream >> trash; //titleY
                linestream >> nbinY;
                for(int i = 0; i <= nbinY; ++i){
                    Double_t temp;
                    linestream>>temp;
                    y_binsArr.push_back(temp);
                }
                linestream >> trash; //#
                linestream >> trash; //titleY
                linestream >> nbinX;
                for(int i = 0; i <= nbinX; ++i){
                    Double_t temp;
                    linestream>>temp;
                    x_binsArr.push_back(temp);
                }
                linestream >> trash; //#
            }
        }
        myfile.close();
     }
     
     ngBins=(nbinX+2)*(nbinY+2);
}





AnalyzerDoubleDiffXS::AnalyzerDoubleDiffXS(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("dda_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerDoubleDiffXS::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    
    //FIXME: Why this function is executing tow times?
    
    TString name;
    int nbinY;
    int nbinX;
    std::vector<Double_t> y_binsArr;
    std::vector<Double_t> x_binsArr;
    int ngBins=0;
    
    //top_rapidity_vs_top_pt
     name = "top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"Top Rapidity vs Top pT; p_{T}^{t} [GeV];y(t)",1200,0,1200,200,-5,5));
     setBinning(name,nbinX,nbinY, x_binsArr, y_binsArr,ngBins); 
     name = "gen_top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"Top Rapidity vs Top pT; p_{T}^{t} [GeV];y(t)",nbinX,&(x_binsArr.at(0)),nbinY,&(y_binsArr.at(0)) ));
     name = "gen_if_reco_top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"Top Rapidity vs Top pT; p_{T}^{t} [GeV];y(t)",nbinX,&(x_binsArr.at(0)),nbinY,&(y_binsArr.at(0)) ));
     name = "genreco_map_top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"global reco vs gen;reco gBins;gen gBins",ngBins,-0.5,ngBins-0.5,ngBins,-0.5,ngBins-0.5 ));
     
    //ttbar_mass_vs_ttbar_rapidity
    name = "ttbar_mass_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs ttbar rapidity;y(ttbar);mass(ttbar), [GeV]",200,-5,5,2000,0,2000));
    setBinning(name,nbinX,nbinY, x_binsArr, y_binsArr,ngBins); 
    name = "gen_ttbar_mass_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs ttbar rapidity;y(ttbar);mass(ttbar), [GeV]",nbinX,&(x_binsArr.at(0)),nbinY,&(y_binsArr.at(0)) ));
    name = "gen_if_reco_ttbar_mass_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs ttbar rapidity;y(ttbar);mass(ttbar), [GeV]",nbinX,&(x_binsArr.at(0)),nbinY,&(y_binsArr.at(0)) ));
    name = "genreco_map_ttbar_mass_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D (prefix_+name+step,"global reco vs gen;reco gBins;gen gBins",ngBins,-0.5,ngBins-0.5,ngBins,-0.5,ngBins-0.5 ));
     
     
     name = "ttbar_rapidity_vs_ttbar_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"TTBar Rapidity vs TTBar pT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",1200,0,1200,200,-5,5));
     name = "gen_ttbar_rapidity_vs_ttbar_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"TTBar Rapidity vs TTBar pT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",1200,0,1200,200,-5,5));
     
    name = "bjet_multiplicity_vs_x";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"B-Jet Multiplicity vs x;;N b-jets",1000,0,1,21,-0.5,20.5));
    name = "gen_bjet_multiplicity_vs_x";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"B-Jet Multiplicity vs x;;N b-jets",1000,0,1,21,-0.5,20.5));
    
    name = "jet_multiplicity_vs_x";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs x;;N b-jets",1000,0,1,21,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_x";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs x;;N b-jets",1000,0,1,21,-0.5,20.5));
    
    name = "dummy_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"dummy vs x1;;",400,0,1,1,0,10));
    name = "dummy_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"dummy vs x2;;",400,0,1,1,0,10));
    
    name = "jet_multiplicity_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top pt;pt(top), [GeV];N b-jets",1200,0,1200,21000,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top pt;pt(top), [GeV];N b-jets",1200,0,1200,21000,-0.5,20.5));
    
    name = "jet_multiplicity_vs_top_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top y;y(top) ;N b-jets",200,-5,5,21,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_top_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top y;y(top) ;N b-jets",200,-5,5,21,-0.5,20.5));
    
    name = "jet_multiplicity_vs_ttbar_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar pt;pt(ttbar), [GeV];N b-jets",1200,0,1200,21,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_ttbar_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar pt;pt(ttbar), [GeV];N b-jets",1200,0,1200,21,-0.5,20.5));
    
    name = "jet_multiplicity_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar y;y(ttbar) ;N b-jets",200,-5,5,21,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_ttbar_rapidity";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar y;y(ttbar) ;N b-jets",200,-5,5,21,-0.5,20.5));
    
    name = "jet_multiplicity_vs_ttbar_mass";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar mass;mass(ttbar), [GeV];N b-jets",2000,0,2000,21,-0.5,20.5));
    name = "gen_jet_multiplicity_vs_ttbar_mass";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar mass;mass(ttbar), [GeV];N b-jets",2000,0,2000,21,-0.5,20.5));
    
    name = "ttbar_mass_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs top pt;pt(top), [GeV];mass(ttbar), [GeV]",1200,0,1200,2000,0,2000));
    name = "gen_ttbar_mass_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs top pt;pt(top), [GeV];mass(ttbar), [GeV]",1200,0,1200,2000,0,2000));
    
    name = "ttbar_delta_phi_vs_ttbar_pt";
    m_histogram[name] = this->store(new TH2D (prefix_+name+step,"#Delta#phi vs TTBar pT;p_{T}^{t#bar{t}} [GeV];#delta(#phi(t#bar{t}))",1200,0,1200,1000,0,3.15));
    name = "gen_ttbar_delta_phi_vs_ttbar_pt";
    m_histogram[name] = this->store(new TH2D (prefix_+name+step,"#Delta#phi vs TTBar pT;p_{T}^{t#bar{t}} [GeV];#delta(#phi(t#bar{t}))",1200,0,1200,1000,0,3.15));
    
    name = "ttbar_delta_phi_vs_ttbar_mass";
    m_histogram[name] = this->store(new TH2D (prefix_+name+step,"#Delta#phi vs ttbar mass;mass(ttbar), [GeV];#delta(#phi(t#bar{t}))",2000,0,2000,1000,0,3.15));
    name = "gen_ttbar_delta_phi_vs_ttbar_mass";
    m_histogram[name] = this->store(new TH2D (prefix_+name+step,"#Delta#phi vs ttbar mass;mass(ttbar), [GeV];#delta(#phi(t#bar{t}))",2000,0,2000,1000,0,3.15));
    
    name = "ttbar_delta_eta_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"delta eta vs top pt;pt(top), [GeV];#delta(#eta(t#bar{t}))",1200,0,1200,1000,0,10));
    name = "gen_ttbar_delta_eta_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"delta eta vs top pt;pt(top), [GeV];#delta(#eta(t#bar{t}))",1200,0,1200,1000,0,10));
    
    
    
//     name = "ttbar_mass_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name,"ttbar mass vs x;x;N b-jets",400,0,1,2000,0,2000));
//     name = "ttbar_pt_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,400,0,400));
//     name = "ttbar_rapidity_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,200,-5,5));
//     name = "top_pt_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,1200,0,1200));
//     name = "top_rapidity_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,200,-5,5));
    
}



void AnalyzerDoubleDiffXS::fillHistos(const EventMetadata&,
                                      const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                      const TopGenObjects& topGenObjects,
                                      const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{

    
   if(kinematicReconstructionSolutions.numberOfSolutions()){
    
        const LV& hypTop = kinematicReconstructionSolutions.solution().top();
        const LV& hypAntiTop = kinematicReconstructionSolutions.solution().antiTop();
        const LV& hypTtbar = kinematicReconstructionSolutions.solution().ttbar();

   //proton Energy [GeV]  
   double protonE = 4000;
   TLorentzVector hyptop(common::LVtoTLV(hypTop));
   TLorentzVector hypantitop(common::LVtoTLV(hypAntiTop));
   TLorentzVector hypttbar(hyptop+hypantitop);
   
   int nRecoJets=recoObjectIndices.jetIndices_.size();
   
   ((TH2D*)m_histogram["top_rapidity_vs_top_pt"])->Fill(hyptop.Pt(),hyptop.Rapidity(),weight);
   //((TH2D*)m_histogram["top_rapidity_vs_top_pt"])->Fill(hypantitop.Pt(),hypantitop.Rapidity(),weight);
   
   ((TH2D*)m_histogram["ttbar_rapidity_vs_ttbar_pt"])->Fill(hypttbar.Pt(),hypttbar.Rapidity(),weight);
   
   double restPzJetsSum=0; // rest means jets not from top or anti-top
   double restEJetsSum=0; 
   for(int i=0;i<nRecoJets;i++)
   {
       if(i != kinematicReconstructionSolutions.solution().bjetIndex() && i != kinematicReconstructionSolutions.solution().bjetIndex())
       {
           restEJetsSum += (*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(i)).E();
           restPzJetsSum += (*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(i)).Pz();
       }
   }
   
   
   double x1 = (hyptop.E()+hypantitop.E()+restEJetsSum+hyptop.Pz()+hypantitop.Pz()+restPzJetsSum)/(2*protonE);
   double x2 = (hyptop.E()+hypantitop.E()+restEJetsSum-hyptop.Pz()-hypantitop.Pz()-restPzJetsSum)/(2*protonE);
   ((TH2D*)m_histogram["bjet_multiplicity_vs_x"])->Fill(x1,recoObjectIndices.bjetIndices_.size(),weight);
   ((TH2D*)m_histogram["bjet_multiplicity_vs_x"])->Fill(x2,recoObjectIndices.bjetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_x"])->Fill(x1,nRecoJets,weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_x"])->Fill(x2,nRecoJets,weight);
   
   
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_pt"])->Fill(hyptop.Pt(),nRecoJets,weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_pt"])->Fill(hypantitop.Pt(),nRecoJets,weight);
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_rapidity"])->Fill(hyptop.Rapidity(),nRecoJets,weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_rapidity"])->Fill(hypantitop.Rapidity(),nRecoJets,weight);
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_pt"])->Fill((hyptop+hypantitop).Pt(),nRecoJets,weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_rapidity"])->Fill((hyptop+hypantitop).Rapidity(),nRecoJets,weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_mass"])->Fill((hyptop+hypantitop).M(),nRecoJets,weight);
   
   ((TH2D*)m_histogram["ttbar_mass_vs_top_pt"])->Fill(hyptop.Pt(),(hyptop+hypantitop).M(),weight);
   ((TH2D*)m_histogram["ttbar_mass_vs_top_pt"])->Fill(hypantitop.Pt(),(hyptop+hypantitop).M(),weight);
   
   ((TH2D*)m_histogram["ttbar_mass_vs_ttbar_rapidity"])->Fill((hyptop+hypantitop).Rapidity(),(hyptop+hypantitop).M(),weight);
   
   ((TH2D*)m_histogram["ttbar_delta_eta_vs_top_pt"])->Fill(hyptop.Pt(),fabs(hyptop.Eta()-hypantitop.Eta()),weight);
   ((TH2D*)m_histogram["ttbar_delta_eta_vs_top_pt"])->Fill(hypantitop.Pt(),fabs(hyptop.Eta()-hypantitop.Eta()),weight);
   
   ((TH2D*)m_histogram["ttbar_delta_phi_vs_ttbar_pt"])->Fill(hypttbar.Pt(),fabs(hyptop.DeltaPhi(hypantitop)),weight);
   
   ((TH2D*)m_histogram["ttbar_delta_phi_vs_ttbar_mass"])->Fill(hypttbar.M(),fabs(hyptop.DeltaPhi(hypantitop)),weight);

   
   }
   
   if(topGenObjects.valuesSet_){
    double trueLevelWeight = genLevelWeights.trueLevelWeight_;
    
    TLorentzVector gentop(common::LVtoTLV((*topGenObjects.GenTop_)));
    TLorentzVector genantitop(common::LVtoTLV((*topGenObjects.GenAntiTop_)));
    TLorentzVector genttbar(gentop+genantitop);
    
    
    ((TH2D*)m_histogram["gen_top_rapidity_vs_top_pt"])->Fill(gentop.Pt(),gentop.Rapidity(),trueLevelWeight);
    //((TH2D*)m_histogram["gen_top_rapidity_vs_top_pt"])->Fill(genantitop.Pt(),genantitop.Rapidity(),trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_ttbar_rapidity_vs_ttbar_pt"])->Fill(genttbar.Pt(),genttbar.Rapidity(),trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_ttbar_mass_vs_top_pt"])->Fill(gentop.Pt(),genttbar.M(),trueLevelWeight);
    ((TH2D*)m_histogram["gen_ttbar_mass_vs_top_pt"])->Fill(genantitop.Pt(),genttbar.M(),trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_ttbar_delta_eta_vs_top_pt"])->Fill(gentop.Pt(),fabs(gentop.Eta()-genantitop.Eta()),trueLevelWeight);
    ((TH2D*)m_histogram["gen_ttbar_delta_eta_vs_top_pt"])->Fill(genantitop.Pt(),fabs(gentop.Eta()-genantitop.Eta()),trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_ttbar_delta_phi_vs_ttbar_pt"])->Fill(genttbar.Pt(),fabs(gentop.DeltaPhi(genantitop)),trueLevelWeight);
    ((TH2D*)m_histogram["gen_ttbar_delta_phi_vs_ttbar_mass"])->Fill(genttbar.M(),fabs(gentop.DeltaPhi(genantitop)),trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_ttbar_mass_vs_ttbar_rapidity"])->Fill((gentop+genantitop).Rapidity(),(gentop+genantitop).M(),weight);
    
    int nGenVisJets = genObjectIndices.genVisJetIndices_.size();
    
    //((TH2D*)m_histogram["gen_jet_multiplicity_vs_x"]->Fill(,nGenVisJets,trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_top_pt"])->Fill(gentop.Pt(),nGenVisJets,trueLevelWeight);
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_top_pt"])->Fill(genantitop.Pt(),nGenVisJets,trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_top_rapidity"])->Fill(gentop.Rapidity(),nGenVisJets,trueLevelWeight);
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_top_rapidity"])->Fill(genantitop.Rapidity(),nGenVisJets,trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_ttbar_pt"])->Fill(genttbar.Pt(),nGenVisJets,trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_ttbar_rapidity"])->Fill(genttbar.Rapidity(),nGenVisJets,trueLevelWeight);
    
    ((TH2D*)m_histogram["gen_jet_multiplicity_vs_ttbar_mass"])->Fill(genttbar.M(),nGenVisJets,trueLevelWeight);
    
   }
   
   if(topGenObjects.valuesSet_ && kinematicReconstructionSolutions.numberOfSolutions()){
       
       TLorentzVector hyptop(common::LVtoTLV(kinematicReconstructionSolutions.solution().top()));
       TLorentzVector hypantitop(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop()));
       TLorentzVector hypttbar(hyptop+hypantitop);
       
       TLorentzVector gentop(common::LVtoTLV((*topGenObjects.GenTop_)));
       TLorentzVector genantitop(common::LVtoTLV((*topGenObjects.GenAntiTop_)));
       TLorentzVector genttbar(gentop+genantitop);
       
       
       int gBinReco=0;
       int gBinGen=0;
       
       
       //top_rapidity_vs_top_pt
       ((TH2D*)m_histogram["gen_if_reco_top_rapidity_vs_top_pt"])->Fill(gentop.Pt(),gentop.Rapidity(),weight);
       //((TH2D*)m_histogram["gen_if_reco_top_rapidity_vs_top_pt"])->Fill(genantitop.Pt(),genantitop.Rapidity(),weight);
       
       gBinReco=m_histogram["gen_top_rapidity_vs_top_pt"]->FindBin(hyptop.Pt(),hyptop.Rapidity());
       gBinGen=m_histogram["gen_top_rapidity_vs_top_pt"]->FindBin(gentop.Pt(),gentop.Rapidity());
       ((TH2D*)m_histogram["genreco_map_top_rapidity_vs_top_pt"])->Fill(gBinReco,gBinGen,weight);
       gBinReco=m_histogram["gen_top_rapidity_vs_top_pt"]->FindBin(hypantitop.Pt(),hypantitop.Rapidity());
       gBinGen=m_histogram["gen_top_rapidity_vs_top_pt"]->FindBin(genantitop.Pt(),genantitop.Rapidity());
       //((TH2D*)m_histogram["genreco_map_top_rapidity_vs_top_pt"])->Fill(gBinReco,gBinGen,weight);
       
       
       //ttbar_mass_vs_ttbar_rapidity
       ((TH2D*)m_histogram["gen_if_reco_ttbar_mass_vs_ttbar_rapidity"])->Fill(genttbar.M(),genttbar.Rapidity(),weight);
       
       gBinReco=m_histogram["gen_ttbar_mass_vs_ttbar_rapidity"]->FindBin(hypttbar.M(),hypttbar.Rapidity());
       gBinGen=m_histogram["gen_ttbar_mass_vs_ttbar_rapidity"]->FindBin(genttbar.M(),genttbar.Rapidity());
       ((TH2D*)m_histogram["genreco_map_ttbar_mass_vs_ttbar_rapidity"])->Fill(gBinReco,gBinGen,weight);
       
//        
   }
   
}








