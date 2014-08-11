#include <iostream>
#include <cstdlib>

#include <Math/VectorUtil.h>
#include "TLorentzVector.h"

#include "VariablesTTBar.h"
#include "analysisStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"





// ----------------------------------------- Methods for VariablesTTBar -----------------------------------------------------------



VariablesTTBar::VariablesTTBar():
VariablesBase(),
top_pt_(VariableFloat(name_top_pt_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_pt_(VariableFloat(name_ttbar_pt_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),
jet_multiplicity_(VariableInt(name_jet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),
gen_jet_multiplicity_(VariableInt(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),

isTopGen_(VariableInt(name_isTopGen_)),
isKinReco_(VariableInt(name_isKinReco_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))

{}



VariablesTTBar::VariablesTTBar(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinRecoObjects& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight):
VariablesBase(weight),
top_pt_(VariableFloat(name_top_pt_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_pt_(VariableFloat(name_ttbar_pt_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),
jet_multiplicity_(VariableInt(name_jet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),
gen_jet_multiplicity_(VariableInt(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),

isTopGen_(VariableInt(name_isTopGen_)),
isKinReco_(VariableInt(name_isKinReco_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))
{
    
    isKinReco_.value_ = 0;
    isTopGen_.value_ = 0;
    
  if(kinRecoObjects.valuesSet_){
   
   isKinReco_.value_ = 1;

   
   //proton Energy [GeV]  
   double protonE = 4000;
   TLorentzVector hyptop(common::LVtoTLV((*kinRecoObjects.HypTop_).at(0)));
   TLorentzVector hypantitop(common::LVtoTLV((*kinRecoObjects.HypAntiTop_).at(0)));
   TLorentzVector hypttbar(hyptop+hypantitop);
   
   int nRecoJets=recoObjectIndices.jetIndices_.size();
   jet_multiplicity_.value_ = nRecoJets;
   
   top_pt_.value_ = hyptop.Pt();
   top_rapidity_.value_ = hyptop.Rapidity();
   ttbar_pt_.value_ = hypttbar.Pt();
   ttbar_rapidity_.value_ = hypttbar.Rapidity();
   ttbar_mass_.value_ = hypttbar.M();
   
   ttbar_delta_eta_.value_ = fabs(hyptop.Eta()-hypantitop.Eta());
   ttbar_delta_phi_.value_ = fabs(hyptop.Eta()-hypantitop.Eta());
   
   double restPzJetsSum=0; // rest means jets not from top or anti-top
   double restEJetsSum=0; 
   for(int i=0;i<nRecoJets;i++)
   {
       if(i != (*kinRecoObjects.HypJet0index_).at(0) && i != (*kinRecoObjects.HypJet1index_).at(0))
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
       
      isTopGen_.value_ = 1;
      trueLevelWeight_.value_ = genLevelWeights.trueLevelWeight_;
      
      gen_jet_multiplicity_.value_ = genObjectIndices.genVisJetIndices_.size();
    
      TLorentzVector gentop(common::LVtoTLV((*topGenObjects.GenTop_)));
      TLorentzVector genantitop(common::LVtoTLV((*topGenObjects.GenAntiTop_)));
      TLorentzVector genttbar(gentop+genantitop);
    
      gen_top_pt_.value_ = gentop.Pt();
      gen_top_rapidity_.value_ = gentop.Rapidity();
      gen_ttbar_pt_.value_ = genttbar.Pt();
      gen_ttbar_rapidity_.value_ = genttbar.Rapidity();
      gen_ttbar_mass_.value_ = genttbar.M();
    
      gen_ttbar_delta_eta_.value_ =fabs(gentop.Eta()-genantitop.Eta());
      gen_ttbar_delta_phi_.value_ =fabs(gentop.DeltaPhi(genantitop));
    
   }
    
}



std::vector<VariablesBase*> VariablesTTBar::fillVariables(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinRecoObjects& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight)
{
    std::vector<VariablesBase*> result;

        VariablesBase* variables = new VariablesTTBar(recoObjects, commonGenObjects, topGenObjects, kinRecoObjects, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights,weight);
        result.push_back(variables);
    
    return result;
}







