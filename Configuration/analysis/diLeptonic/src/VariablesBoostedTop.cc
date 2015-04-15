#include <iostream>
#include <cstdlib>

#include <Math/VectorUtil.h>
#include "TLorentzVector.h"

#include "VariablesBoostedTop.h"
#include "analysisStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"




// ----------------------------------------- Methods for VariablesBoostedTop -----------------------------------------------------------



VariablesBoostedTop::VariablesBoostedTop():
VariablesBase(),
lep_pt_(VariableFloat(name_lep_pt_)),
anti_lep_pt_(VariableFloat(name_anti_lep_pt_)),

top_pt_(VariableFloat(name_top_pt_)),
topbar_pt_(VariableFloat(name_topbar_pt_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_pt_(VariableFloat(name_ttbar_pt_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),
jet_multiplicity_(VariableInt(name_jet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),
mlblbmet_(VariableFloat(name_mlblbmet_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_topbar_pt_(VariableFloat(name_gen_topbar_pt_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),
gen_jet_multiplicity_(VariableInt(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),
gen_mlblbmet_(VariableFloat(name_gen_mlblbmet_)),

entry_(VariableInt(name_entry_)),
isTopGen_(VariableInt(name_isTopGen_)),
isKinReco_(VariableInt(name_isKinReco_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))

{}



VariablesBoostedTop::VariablesBoostedTop(const EventMetadata& eventMetadata,
                               const RecoObjects& recoObjects, const CommonGenObjects& ,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& ,
                                                          const double& weight):
VariablesBase(weight),
lep_pt_(VariableFloat(name_lep_pt_)),
anti_lep_pt_(VariableFloat(name_anti_lep_pt_)),

top_pt_(VariableFloat(name_top_pt_)),
topbar_pt_(VariableFloat(name_topbar_pt_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_pt_(VariableFloat(name_ttbar_pt_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),
jet_multiplicity_(VariableInt(name_jet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),
mlblbmet_(VariableFloat(name_mlblbmet_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_topbar_pt_(VariableFloat(name_gen_topbar_pt_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),
gen_jet_multiplicity_(VariableInt(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),
gen_mlblbmet_(VariableFloat(name_gen_mlblbmet_)),

entry_(VariableInt(name_entry_)),
isTopGen_(VariableInt(name_isTopGen_)),
isKinReco_(VariableInt(name_isKinReco_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))
{
    
    isKinReco_.value_ = 0;
    isTopGen_.value_ = 0;
    entry_.value_ = -999;
    

    
  if(kinematicReconstructionSolutions.numberOfSolutions()){
   
   isKinReco_.value_ = 1;
    
   //proton Energy [GeV]  
   double protonE = 4000;
   TLorentzVector hyptop(common::LVtoTLV(kinematicReconstructionSolutions.solution().top()));
   TLorentzVector hypantitop(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop()));
   TLorentzVector hypbjet(common::LVtoTLV(kinematicReconstructionSolutions.solution().bjet()));
   TLorentzVector hypantibjet(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiBjet()));
   TLorentzVector hyplep(common::LVtoTLV(kinematicReconstructionSolutions.solution().lepton()));
   TLorentzVector hypantilep(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiLepton()));
   TLorentzVector hypn(common::LVtoTLV(kinematicReconstructionSolutions.solution().neutrino()));
   TLorentzVector hypantin(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiNeutrino()));
   TLorentzVector hypttbar(hyptop+hypantitop);
   TLorentzVector hypmet((hypn+hypantin).X(),(hypn+hypantin).Y(),0,0);
   
   int nRecoJets=recoObjectIndices.jetIndices_.size();
   jet_multiplicity_.value_ = nRecoJets;
   
   lep_pt_.value_      = hyplep.Pt();
   anti_lep_pt_.value_ = hypantilep.Pt();
   
   top_pt_.value_ = hyptop.Pt();
   topbar_pt_.value_ = hypantitop.Pt();
   top_rapidity_.value_ = hyptop.Rapidity();
   ttbar_pt_.value_ = hypttbar.Pt();
   ttbar_rapidity_.value_ = hypttbar.Rapidity();
   ttbar_mass_.value_ = hypttbar.M();
   
   ttbar_delta_eta_.value_ = fabs(hyptop.Eta()-hypantitop.Eta());
   ttbar_delta_phi_.value_ = fabs(hyptop.Eta()-hypantitop.Eta());
   
   mlblbmet_.value_ = (hypmet+hyplep+hypantilep+hypbjet+hypantibjet).M()/hypttbar.M();
   
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
   x1_.value_ = x1;
   x2_.value_ = x2;
   
   }
   
   if(topGenObjects.valuesSet_){
       
      entry_.value_ = (Int_t)eventMetadata.eventNumber_;
      isTopGen_.value_ = 1;
      trueLevelWeight_.value_ = genLevelWeights.trueLevelWeight_;
      
      gen_jet_multiplicity_.value_ = genObjectIndices.genVisJetIndices_.size();
    
      TLorentzVector gentop(common::LVtoTLV((*topGenObjects.GenTop_)));
      TLorentzVector genantitop(common::LVtoTLV((*topGenObjects.GenAntiTop_)));
      TLorentzVector genbjet(common::LVtoTLV((*topGenObjects.GenB_)));
      TLorentzVector genantibjet(common::LVtoTLV((*topGenObjects.GenAntiB_)));
      TLorentzVector genlep(common::LVtoTLV((*topGenObjects.GenLepton_)));
      TLorentzVector genantilep(common::LVtoTLV((*topGenObjects.GenAntiLepton_)));
      TLorentzVector genn(common::LVtoTLV((*topGenObjects.GenNeutrino_)));
      TLorentzVector genantin(common::LVtoTLV((*topGenObjects.GenAntiNeutrino_)));
      TLorentzVector genttbar(gentop+genantitop);
      TLorentzVector genmet((genn+genantin).X(),(genn+genantin).Y(),0,0);
    
      gen_top_pt_.value_ = gentop.Pt();
      gen_topbar_pt_.value_ = genantitop.Pt();
      gen_top_rapidity_.value_ = gentop.Rapidity();
      gen_ttbar_pt_.value_ = genttbar.Pt();
      gen_ttbar_rapidity_.value_ = genttbar.Rapidity();
      gen_ttbar_mass_.value_ = genttbar.M();
    
      gen_ttbar_delta_eta_.value_ =fabs(gentop.Eta()-genantitop.Eta());
      gen_ttbar_delta_phi_.value_ =fabs(gentop.DeltaPhi(genantitop));
    
      gen_mlblbmet_.value_ = (genmet+genlep+genantilep+genbjet+genantibjet).M()/genttbar.M();
      
   }
    
}



std::vector<VariablesBase*> VariablesBoostedTop::fillVariables(const EventMetadata& eventMetadata, 
                                                          const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight)
{
    std::vector<VariablesBase*> result;

        VariablesBase* variables = new VariablesBoostedTop(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights,weight);
        result.push_back(variables);
    
    return result;
}







