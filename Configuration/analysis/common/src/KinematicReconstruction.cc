#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include <TLorentzVector.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <Math/VectorUtil.h>

#include "classes.h"
#include "utils.h"
#include "analysisUtils.h"
#include "KinematicReconstruction.h"
#include "KinematicReconstruction_LSroutines.h"
#include "KinematicReconstruction_MeanSol.h"
#include "KinematicReconstructionSolution.h"




constexpr double TopMASS = 172.5;













// -------------------------------------- Methods for KinematicReconstruction --------------------------------------

void KinematicReconstruction::angle_rot(const double& alpha, const double& e, const TLorentzVector& inJet, TLorentzVector& jet_sm)const
{
    double px_1, py_1, pz_1; // Coordinate system where momentum is along Z-axis
    
    //Transition matrix detector -> syst1 ...
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    // ...
    
    TLorentzVector jet = inJet;
    if(fabs(jet.Px())<=e){jet.SetPx(0);}
    if(fabs(jet.Py())<=e){jet.SetPy(0);}
    if(fabs(jet.Pz())<=e){jet.SetPz(0);}
    
    //Rotation in syst 1 ...
    double phi = 2*TMath::Pi()*r3_->Rndm();
    pz_1 = jet.Vect().Mag()*cos(alpha);
    px_1 = - jet.Vect().Mag()*sin(alpha)*sin(phi);
    py_1 = jet.Vect().Mag()*sin(alpha)*cos(phi);  
    // ...
    
    //Transition detector <- syst1 ...
    if (jet.Py()!=0||jet.Pz()!=0)
    {
        double d = sqrt(jet.Pz()*jet.Pz() + jet.Py()*jet.Py());
        double p = jet.Vect().Mag();
        
        x1 = d/p;
        y1 = 0;
        z1 = jet.Px()/p;
        
        x2 = - jet.Px()*jet.Py()/d/p;
        y2 = jet.Pz()/d;
        z2 = jet.Py()/p;
        
        x3 = - jet.Px()*jet.Pz()/d/p;
        y3 = - jet.Py()/d;
        z3 = jet.Pz()/p;
        
        jet_sm.SetPx(x1*px_1+y1*py_1+z1*pz_1);
        jet_sm.SetPy(x2*px_1+y2*py_1+z2*pz_1);
        jet_sm.SetPz(x3*px_1+y3*py_1+z3*pz_1);
        jet_sm.SetE(jet.E());
    }
    
    if (jet.Px()==0&&jet.Py()==0&&jet.Pz()==0)
    {
        jet_sm.SetPx(jet.Px());
        jet_sm.SetPy(jet.Py());
        jet_sm.SetPz(jet.Pz());
        jet_sm.SetE(jet.E());
    }
    
    if (jet.Px()!=0&&jet.Py()==0&&jet.Pz()==0)
    {
        jet_sm.SetPx(pz_1);
        jet_sm.SetPy(px_1);
        jet_sm.SetPz(py_1);
        jet_sm.SetE(jet.E());
    }
    // ...
}



KinematicReconstruction::KinematicReconstruction(const int minNumberOfBtags, const bool preferBtags, const bool massLoop):
minNumberOfBtags_(minNumberOfBtags),
preferBtags_(preferBtags),
massLoop_(massLoop),
nSol_(0),
isJetsMerging_(false),
h_wmass_(0),
h_jetAngleRes_(0),
h_jetEres_(0),
h_lepAngleRes_(0),
h_lepEres_(0),
h_mbl_w_(0)
{
    std::cout<<"--- Beginning preparation of kinematic reconstruction\n";
    
    // Sanity check
    if(minNumberOfBtags_<0 || minNumberOfBtags_>2){
        std::cerr<<"ERROR in constructor of KinematicReconstruction! Minimum number of b-tags needs to be within [0, 2]. Value is: "
                 <<minNumberOfBtags_<<"\n...break\n"<<std::endl;
        exit(96);
    }
    
    std::cout<<"Require minimum number of b-tags per solution: "<<minNumberOfBtags_<<"\n";
    std::cout<<"Prefer solutions with more b-tags: "<<(preferBtags_ ? "yes" : "no")<<"\n";
    std::cout<<"Uncertainty treatment: "<<(massLoop_ ? "top mass scan" : "random-number based smearing")<<"\n";
    
    // Read all histograms for smearings from file
    if(!massLoop_) this->loadData();
    
    std::cout<<"=== Finishing preparation of kinematic reconstruction\n\n";
}



void KinematicReconstruction::setRandomNumberSeeds(const LV& lepton, const LV& antiLepton, const LV& jet1, const LV& jet2)const
{
    // Asymmetric treatment of both jets and also both leptons, to ensure different seed for each combination
    const unsigned int seed =  std::abs(static_cast<int>( 1.e6*(jet1.pt()/jet2.pt()) * std::sin((lepton.Pt() + 2.*antiLepton.pt())*1.e6) ) );
    gRandom->SetSeed(seed);
    r3_->SetSeed(seed);
}



KinematicReconstructionSolutions KinematicReconstruction::solutions(const std::vector<int>& leptonIndices, const std::vector<int>& antiLeptonIndices,
                                                                    const std::vector<int>& jetIndices, const std::vector<int>& bjetIndices,
                                                                    const VLV& allLeptons,
                                                                    const VLV& allJets, const std::vector<double>& btags,
                                                                    const LV& met)const
{
    KinematicReconstructionSolutions result;
    
    // Check if minimum number of objects for any solution is present
    const int numberOfJets(jetIndices.size());
    const int numberOfBtags(bjetIndices.size());
    if(!leptonIndices.size() || !antiLeptonIndices.size() || numberOfJets<2 || numberOfBtags<minNumberOfBtags_) return result;
    
    // Find solutions with 2 b-tagged jets
    for(std::vector<int>::const_iterator i_index = bjetIndices.begin(); i_index != bjetIndices.end(); ++i_index)
        for(std::vector<int>::const_iterator j_index = i_index+1; j_index != bjetIndices.end(); ++j_index)
            for(const int leptonIndex : leptonIndices)
                for(const int antiLeptonIndex : antiLeptonIndices){
                    const std::vector<KinematicReconstructionSolution> solutions(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, *i_index, *j_index, allLeptons, allJets, btags, met, 2));
                    result.addSolutions(solutions);
                    const std::vector<KinematicReconstructionSolution> solutionsSwapped(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, *j_index, *i_index, allLeptons, allJets, btags, met, 2));
                    result.addSolutions(solutionsSwapped);
                }
    if(preferBtags_ && result.numberOfSolutionsTwoBtags()) return result;
    if(minNumberOfBtags_ > 1) return result;
    
    // Indices of non b-tagged jets
    const int numberOfNonBtags(numberOfJets - numberOfBtags);
    if(numberOfNonBtags < 1) return result;
    std::vector<int> nonBjetIndices;
    for(const int index : jetIndices)
        if(std::find(bjetIndices.begin(), bjetIndices.end(), index) == bjetIndices.end()) nonBjetIndices.push_back(index);
    
    // Find solutions with 1 b-tagged jet
    for(const int bjetIndex : bjetIndices)
        for(const int nonBjetIndex : nonBjetIndices)
            for(const int leptonIndex : leptonIndices)
                for(const int antiLeptonIndex : antiLeptonIndices){
                    const std::vector<KinematicReconstructionSolution> solutions(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, bjetIndex, nonBjetIndex, allLeptons, allJets, btags, met, 1));
                    result.addSolutions(solutions);
                    const std::vector<KinematicReconstructionSolution> solutionsSwapped(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, nonBjetIndex, bjetIndex, allLeptons, allJets, btags, met, 1));
                    result.addSolutions(solutionsSwapped);
                }
    if(preferBtags_ && result.numberOfSolutionsOneBtag()) return result;
    if(minNumberOfBtags_ > 0) return result;
    
    // Find solutions with 0 b-tagged jets
    for(std::vector<int>::const_iterator i_index = nonBjetIndices.begin(); i_index != nonBjetIndices.end(); ++i_index)
        for(std::vector<int>::const_iterator j_index = i_index+1; j_index != nonBjetIndices.end(); ++j_index)
            for(const int leptonIndex : leptonIndices)
                for(const int antiLeptonIndex : antiLeptonIndices){
                    const std::vector<KinematicReconstructionSolution> solutions(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, *i_index, *j_index, allLeptons, allJets, btags, met, 0));
                    result.addSolutions(solutions);
                    const std::vector<KinematicReconstructionSolution> solutionsSwapped(this->solutionsPerObjectCombination(leptonIndex, antiLeptonIndex, *j_index, *i_index, allLeptons, allJets, btags, met, 0));
                    result.addSolutions(solutionsSwapped);
                }
    
    return result;
}



std::vector<KinematicReconstructionSolution> KinematicReconstruction::solutionsPerObjectCombination(const int leptonIndex, const int antiLeptonIndex,
                                                                                                    const int jetIndex1, const int jetIndex2,
                                                                                                    const VLV& allLeptons,
                                                                                                    const VLV& allJets, const std::vector<double>& /*btag values*/ ,
                                                                                                    const LV& met,
                                                                                                    const int numberOfBtags)const
{
    std::vector<KinematicReconstructionSolution> result;
    
    const LV& lepton = allLeptons.at(leptonIndex);
    const LV& antiLepton = allLeptons.at(antiLeptonIndex);
    const LV& jet1 = allJets.at(jetIndex1);
    const LV& jet2 = allJets.at(jetIndex2);
    
    if(massLoop_){
        double neutrinoWeightMax = 0.;
        for(double iTopMass = 100.; iTopMass < 300.5; iTopMass += 1.){
            KinematicReconstruction_LSroutines tp_m(iTopMass, 4.8, 80.4, 0.0, 0.0);
            tp_m.setConstraints(antiLepton, lepton, jet1, jet2, met.px(), met.py());
            if(tp_m.getNsol() < 1) continue;
            if(!(tp_m.getTtSol()->at(0).weight > neutrinoWeightMax)) continue;
            
            neutrinoWeightMax = tp_m.getTtSol()->at(0).weight;
            const LV neutrino = common::TLVtoLV(tp_m.getTtSol()->at(0).neutrino);
            const LV antiNeutrino = common::TLVtoLV(tp_m.getTtSol()->at(0).neutrinobar);
            const LV top = antiLepton + neutrino + jet1;
            const LV antiTop = lepton + antiNeutrino + jet2;
            const double& weight = tp_m.getTtSol()->at(0).weight;
            std::map<KinematicReconstructionSolution::WeightType, double> m_weight;
            m_weight[KinematicReconstructionSolution::defaultForMethod] = weight;
            m_weight[KinematicReconstructionSolution::neutrinoEnergy] = weight;
            const KinematicReconstructionSolution solution(&allLeptons, &allJets,
                                                           leptonIndex, antiLeptonIndex, jetIndex1, jetIndex2,
                                                           top, antiTop, neutrino, antiNeutrino,
                                                           iTopMass, numberOfBtags,
                                                           m_weight);
            result.push_back(solution);
        }
    }
    else{
        KinematicReconstruction_MeanSol meanSolution(TopMASS);
        const bool hasSolution(this->solutionSmearing(meanSolution, lepton, antiLepton, jet1, jet2, met));
        if(hasSolution){
            LV top;
            LV antiTop;
            LV neutrino;
            LV antiNeutrino;
            meanSolution.meanSolution(top, antiTop, neutrino, antiNeutrino);
            const double& weight = meanSolution.getSumWeight();
            std::map<KinematicReconstructionSolution::WeightType, double> m_weight;
            m_weight[KinematicReconstructionSolution::defaultForMethod] = weight;
            m_weight[KinematicReconstructionSolution::averagedSumSmearings_mlb] = weight;
            const KinematicReconstructionSolution solution(&allLeptons, &allJets,
                                                           leptonIndex, antiLeptonIndex, jetIndex1, jetIndex2,
                                                           top, antiTop, neutrino, antiNeutrino,
                                                           TopMASS, numberOfBtags, m_weight);
            result.push_back(solution);
        }
    }
    
    return result;
}



void KinematicReconstruction::doJetsMerging(const VLV* jets, const std::vector<double>* btags)
{
    isJetsMerging_ = true;   
    
    alljets_.clear();
    allbtags_.clear();
    alljets_ = *jets;
    allbtags_ = *btags;
}



void KinematicReconstruction::kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV* jets, const std::vector<double>* btags, const LV* met)
{
    
    sols_.clear();

    //jets selection
    std::vector<int> b1_id;
    std::vector<int> b2_id;
    std::vector<int> nb_tag;
    VLV new_jets(*jets);
    std::vector<double> new_btags(*btags);
    if(!isJetsMerging_) this->inputNoJetMerging(b1_id, b2_id, nb_tag, *btags);
    else this->inputJetMerging(b1_id, b2_id, nb_tag, new_jets, new_btags);
    if(b1_id.size() < 2)return;     

    KinematicReconstruction_MeanSol meanSolution(TopMASS);
    for(int ib = 0; ib < (int)b1_id.size(); ++ib){
        const int bjetIndex = b1_id.at(ib);
        const int antiBjetIndex = b2_id.at(ib);
        const int numberOfBtags = nb_tag.at(ib);
        
        const LV& jet1 = new_jets.at(bjetIndex);
        const LV& jet2 = new_jets.at(antiBjetIndex);
        
        const bool hasSolution(this->solutionSmearing(meanSolution, leptonMinus, leptonPlus, jet1, jet2, *met));
        
        if(hasSolution){
            meanSolution.getMeanSol(sol_.top, sol_.topBar, sol_.neutrino, sol_.neutrinoBar);
            sol_.weight = meanSolution.getSumWeight();
            sol_.Wplus = sol_.lp + sol_.neutrino;
            sol_.Wminus = sol_.lm + sol_.neutrinoBar;
            sol_.ttbar = sol_.top + sol_.topBar;
            sol_.jetB_index = bjetIndex;
            sol_.jetBbar_index = antiBjetIndex;
            sol_.ntags = numberOfBtags;
            sols_.push_back(sol_);
        }
        meanSolution.clear();
    }

    this->setSolutions();
}



bool KinematicReconstruction::solutionSmearing(KinematicReconstruction_MeanSol& meanSolution,
                                               const LV& lepton, const LV& antiLepton,
                                               const LV& jet1, const LV& jet2,
                                               const LV& met)const
{
    // Set random number generator seeds
    this->setRandomNumberSeeds(lepton, antiLepton, jet1, jet2);
    
    const TLorentzVector l_temp = common::LVtoTLV(lepton);
    const TLorentzVector al_temp = common::LVtoTLV(antiLepton);
    const TLorentzVector b_temp = common::LVtoTLV(jet1);
    const TLorentzVector bbar_temp = common::LVtoTLV(jet2);
    const TLorentzVector met_temp = common::LVtoTLV(met);
    
    if((al_temp+b_temp).M()>180. || (l_temp+bbar_temp).M()>180.) return false;
    
    bool isHaveSol(false);
    
//          double b_w=new_btags.at(j1)/(1-new_btags.at(j1));
//          double bbar_w=new_btags.at(j2)/(1-new_btags.at(j2));
//          double b_w=1;
//          double bbar_w=1;
        
        TVector3 vX_reco =  - b_temp.Vect() - bbar_temp.Vect() - l_temp.Vect() - al_temp.Vect() - met_temp.Vect();
        
        // Loop over all smearings
        for(int sm=0; sm<100; ++sm){
            TLorentzVector b_sm=b_temp;
            TLorentzVector bbar_sm=bbar_temp;
            TLorentzVector met_sm;
            TLorentzVector l_sm=l_temp;
            TLorentzVector al_sm=al_temp;
            TLorentzVector vX_sm;
            //jets energy smearing
            double fB=h_jetEres_->GetRandom();//fB=1;  //sm off
            double xB=sqrt((fB*fB*b_sm.E()*b_sm.E()-b_sm.M2())/(b_sm.P()*b_sm.P()));
            double fBbar=h_jetEres_->GetRandom();//fBbar=1; //sm off
            double xBbar=sqrt((fBbar*fBbar*bbar_sm.E()*bbar_sm.E()-bbar_sm.M2())/(bbar_sm.P()*bbar_sm.P()));
            //leptons energy smearing
            double fL=h_lepEres_->GetRandom();//fL=1; //sm off
            double xL=sqrt((fL*fL*l_sm.E()*l_sm.E()-l_sm.M2())/(l_sm.P()*l_sm.P()));
            double faL=h_lepEres_->GetRandom();//faL=1;  //sm off
            double xaL=sqrt((faL*faL*al_sm.E()*al_sm.E()-al_sm.M2())/(al_sm.P()*al_sm.P()));
            //b-jet angle smearing
            b_sm.SetXYZT(b_sm.Px()*xB,b_sm.Py()*xB,b_sm.Pz()*xB,b_sm.E()*fB);
            angle_rot(h_jetAngleRes_->GetRandom(),0.001,b_sm,b_sm);
            //bbar jet angel smearing
            bbar_sm.SetXYZT(bbar_sm.Px()*xBbar,bbar_sm.Py()*xBbar,bbar_sm.Pz()*xBbar,bbar_sm.E()*fBbar);    
            angle_rot(h_jetAngleRes_->GetRandom(),0.001,bbar_sm,bbar_sm);
            //lepton angle smearing
            l_sm.SetXYZT(l_sm.Px()*xL,l_sm.Py()*xL,l_sm.Pz()*xL,l_sm.E()*fL);
            angle_rot(h_lepAngleRes_->GetRandom(),0.001,l_sm,l_sm);
            // anti lepton angle smearing
            al_sm.SetXYZT(al_sm.Px()*xaL,al_sm.Py()*xaL,al_sm.Pz()*xaL,al_sm.E()*faL);
            angle_rot(h_lepAngleRes_->GetRandom(),0.001,al_sm,al_sm);

//X Pt+angle smearing
            double fPx=1;
            double fPy=1;
            double tempX= vX_reco.Px()*fPx;
            double tempY= vX_reco.Py()*fPy;
            vX_sm.SetXYZM(tempX,tempY,0,0);
//...
            
            TVector3 metV3_sm= -b_sm.Vect()-bbar_sm.Vect()-l_sm.Vect()-al_sm.Vect()-vX_sm.Vect();
                met_sm.SetXYZM(metV3_sm.Px(),metV3_sm.Py(),0,0);
            
            KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom());
                tp_sm.setConstraints(al_sm, l_sm, b_sm, bbar_sm, met_sm.Px(), met_sm.Py());

            if(tp_sm.getNsol()>0)
            {
                isHaveSol = true;
                // FIXME: this loop is processed only once by definition, what is it needed for?
                for(int i=0; i<=tp_sm.getNsol()*0; ++i){
                    //meanSolution.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,1);
                    //meanSolution.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,(tp_sm.GetTtSol()->at(i).lepEw));//*b_w*bbar_w
                    
                    double mbl_weight = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((al_sm+b_sm).M()))*h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((l_sm+bbar_sm).M()))/100000000;
                    double pdf_weight = 1;
		   if(0){ // 0 -default , 1 - boosted top
		      double mtt = tp_sm.getTtSol()->at(i).mtt;
		      double x1 = tp_sm.getTtSol()->at(i).x1;
		      double x2 = tp_sm.getTtSol()->at(i).x2;
		      int nMttBins = (int)mttBinsForX_.size()-1;
		      for(int j=0;j<nMttBins;++j)
		      {  
			  if(mtt<mttBinsForX_.at(j)&&j==0){
			      pdf_weight = pdf_weight*h_pdf_w_[0]->GetBinContent(h_pdf_w_[0]->FindBin(x1))*h_pdf_w_[0]->GetBinContent(h_pdf_w_[0]->FindBin(x2));
			  }
			  else if(mtt<mttBinsForX_.at(j+1)&&mtt>mttBinsForX_.at(j))
			  {
			      pdf_weight = pdf_weight*h_pdf_w_[j]->GetBinContent(h_pdf_w_[j]->FindBin(x1))*h_pdf_w_[j]->GetBinContent(h_pdf_w_[j]->FindBin(x2));
			  }
			  else if(mtt>mttBinsForX_.at(j)&&j==nMttBins){
			      pdf_weight = pdf_weight*h_pdf_w_[nMttBins-1]->GetBinContent(h_pdf_w_[nMttBins-1]->FindBin(x1))*h_pdf_w_[nMttBins-1]->GetBinContent(h_pdf_w_[nMttBins-1]->FindBin(x2));
			  } 
		      }
		   }
                    // ...
                    
                    meanSolution.add(tp_sm.getTtSol()->at(i).top,tp_sm.getTtSol()->at(i).topbar,tp_sm.getTtSol()->at(i).neutrino,tp_sm.getTtSol()->at(i).neutrinobar,mbl_weight*pdf_weight);
                    //meanSolution.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,(tp_sm.GetTtSol()->at(i).vw));//
                }
            }
        }

//         if(isHaveSol) {
//             
//                 if(nb_tag[ib]==nbtag){
//                     if(meanSol.getSumWeight()>max_sum_weight) {
//                         max_sum_weight=meanSol.getSumWeight();
//                         meanSol.getMeanSol(sol_.top, sol_.topBar, sol_.neutrino, sol_.neutrinoBar);
//                         sol_.jetB = b_temp;
//                         sol_.jetBbar = bbar_temp;
//                         sol_.lm = l_temp;
//                         sol_.lp = al_temp;
//                         sol_.met = met_temp;
//                         sol_.weight = max_sum_weight;
//                         sol_.Wplus = sol_.lp + sol_.neutrino;
//                         sol_.Wminus = sol_.lm + sol_.neutrinoBar;
//                         sol_.ttbar = sol_.top + sol_.topBar;
//                         sol_.jetB_index = b1_id[ib];
//                         sol_.jetBbar_index = b2_id[ib];
//                         sol_.ntags = nb_tag[ib];
//                         }
//                 }
//             else if(nb_tag[ib]>nbtag) {
//                 nbtag=nb_tag[ib];
//                 max_sum_weight=meanSol.getSumWeight();
//                 meanSol.getMeanSol(sol_.top, sol_.topBar, sol_.neutrino, sol_.neutrinoBar);
//                 sol_.jetB = b_temp;
//                 sol_.jetBbar = bbar_temp;
//                 sol_.lm = l_temp;
//                 sol_.lp = al_temp;
//                 sol_.met = met_temp;
//                 sol_.weight = max_sum_weight;
//                 sol_.Wplus = sol_.lp + sol_.neutrino;
//                 sol_.Wminus = sol_.lm + sol_.neutrinoBar;
//                 sol_.ttbar = sol_.top + sol_.topBar;
//                 sol_.jetB_index = b1_id[ib];
//                 sol_.jetBbar_index = b2_id[ib];
//                 sol_.ntags = nb_tag[ib];
//             }
//         }

    return isHaveSol;
}



void KinematicReconstruction::inputNoJetMerging(std::vector<int>& b1_id, std::vector<int>& b2_id, std::vector<int>& nb_tag,
                                                const std::vector<double>& btags)const
{
    constexpr double btag_wp = 0.244;
    
    for(int i = 0; i < (int)btags.size(); ++i){
        for(int j = 0; j < (int)btags.size(); ++j){
            double wi = btags.at(i);
            double wj = btags.at(j);
            //if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0)) continue;
            if(i==j || (wi<btag_wp && wj<btag_wp)) continue;

            if(wi>btag_wp && wj>btag_wp) nb_tag.push_back(2);
            else nb_tag.push_back(1);

            b1_id.push_back(i);
            b2_id.push_back(j);
        }
    }
}



void KinematicReconstruction::inputJetMerging(std::vector<int>& b1_id, std::vector<int>& b2_id, std::vector<int>& nb_tag,
                                              VLV& new_jets, std::vector<double>& new_btags)const
{
    constexpr double btag_wp = 0.244;
    
    // Fill new jets: logik I.
    //std::vector<TLorentzVector> new_jets_tlv;
    //std::vector<double> new_btag;
    //std::vector<int> veto_i;
    //std::vector<int> veto_j;
    //for(int i=0; i<(int)alljets_.size(); ++i){
    //    if(allbtags_.at(i)<btag_wp)continue;
    //    TLorentzVector ijet=alljets_.at(i);
    //    
    //    for(int j=0; j<(int)alljets_.size(); ++j){
    //        if(i==j)continue;
    //        if(allbtags_.at(j)>=btag_wp)continue; 
    //        TLorentzVector jjet=alljets_.at(j);
    //        double dRij=ijet.DeltaR(jjet);
    //        if(dRij>1)continue;
    //        TLorentzVector new_jet=ijet+jjet;
    //        if(new_jet.Pt()<30)continue;
    //        if(fabs(new_jet.Eta())>2.4)continue;
    //        new_jets_tlv.push_back(new_jet);
    //        new_btag.push_back(allbtags_.at(i));
    //        veto_i.push_back(i);
    //        veto_j.push_back(j);
    //    }   
    //}
    //
    //std::vector<int> index_i;
    //std::vector<int> index_j;
    //for(int i=0; i<(int)alljets_.size(); ++i){
    //    if(alljets_.at(i).Pt()<30)continue;
    //    if(fabs(alljets_.at(i).Eta())>2.4)continue;
    //    jets_tlv.push_back(alljets_.at(i));
    //    new_btags.push_back(allbtags_.at(i));
    //    index_i.push_back(i);
    //    index_j.push_back(i);
    //}
    //
    //for(int i=0; i<(int)new_jets_tlv.size(); ++i){
    //    jets_tlv.push_back(new_jets_tlv.at(i));
    //    new_btags.push_back(new_btag.at(i));
    //    index_i.push_back(veto_i.at(i));
    //    index_j.push_back(veto_j.at(i));
    //}
    //
    //for(int i=0; i<(int)jets_tlv.size(); ++i){
    //    for(int j=0; j<(int)jets_tlv.size(); ++j){
    //        if(index_i[i]==index_i[j] || index_i[i]==index_j[j] || index_j[i]==index_i[j] || index_j[i]==index_j[j])continue;                   
    //        double wi = new_btags.at(i);
    //        double wj = new_btags.at(j);
    //        //if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0))continue;
    //        if(i==j || (wi<btag_wp && wj<btag_wp))continue;
    //        if(wi>btag_wp && wj>btag_wp){nb_tag.push_back(2); }
    //        else{nb_tag.push_back(1); }
    //        b1_id.push_back(i);
    //        b2_id.push_back(j);
    //    }
    //}


    // Fill new jets: logik II.
    new_jets.clear();
    new_btags.clear();
    std::vector<int> index_dR;
    int index_i;
    int index_j;
    
    for(int i = 0; i < (int)alljets_.size(); ++i) index_dR.push_back(-1);
    
    double min_dR(1.);
    do{
        min_dR = 1.;
        index_i = -1;
        index_j = -1;
        for(int i = 0; i < (int)alljets_.size(); ++i){
            if(alljets_.at(i).Pt() < 10) continue;
            if(index_dR.at(i) > -1) continue;
            const LV& ijet = alljets_.at(i);
            
            for(int j = 0; j < (int)alljets_.size(); ++j){
                if(i == j) continue;
                if(alljets_.at(j).Pt() < 10) continue;
                if(index_dR[j]>-1) continue;
                const LV& jjet = alljets_.at(j);
                
                const double dRij = ROOT::Math::VectorUtil::DeltaR(ijet, jjet);
                if(dRij >= 1.) continue;
                if(dRij < min_dR){
                    min_dR = dRij;
                    index_i = i;
                    index_j = j;
                }
            }
        }
        if(index_i >= 0){
            index_dR.at(index_i) = index_j;
            index_dR.at(index_j) = index_i;
            index_i = -1;
            index_j = -1;
        }
    }
    while(min_dR < 1.);
    
    
    std::vector<int> used_index;
    for(int i = 0; i < (int)alljets_.size(); ++i){
        bool flag(false);
        for(int k = 0; k < (int)used_index.size(); ++k) if(used_index.at(k) == i) flag = true;
        if(flag) continue;
        
        if(index_dR[i] < 0){
            if(alljets_.at(i).Pt() < 30.)continue;
            if(fabs(alljets_.at(i).Eta()) > 2.4)continue;
            new_jets.push_back(alljets_.at(i));
            new_btags.push_back(allbtags_.at(i));
        }
        else if(index_dR[i] > -1){
            const LV temp = alljets_.at(i) + alljets_.at((index_dR.at(i)));
            used_index.push_back(index_dR.at(i));
            if(temp.Pt() < 30.) continue;
            if(fabs(temp.Eta()) > 2.4) continue;
            new_jets.push_back(temp);
            new_btags.push_back(allbtags_.at(i)*(allbtags_.at(i)>allbtags_.at(index_dR[i])) + allbtags_.at(index_dR[i])*(allbtags_.at(index_dR[i])>allbtags_.at(i)));
        }
    }
    
    for(int i = 0; i < (int)new_jets.size(); ++i){
        for(int j = 0; j < (int)new_jets.size(); ++j){
            const double wi = new_btags.at(i);
            const double wj = new_btags.at(j);
            //if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0)) continue;
            if(i==j || (wi<btag_wp && wj<btag_wp)) continue;
            if(wi>btag_wp && wj>btag_wp) nb_tag.push_back(2);
            else nb_tag.push_back(1);
            b1_id.push_back(i);
            b2_id.push_back(j);
        }
    }
}



void KinematicReconstruction::setSolutions()
{
    nSol_ = (int)(sols_.size());
    
    if(nSol_ > 0){
        std::nth_element(begin(sols_), begin(sols_), end(sols_),
                         [](const Struct_KinematicReconstruction& a, const Struct_KinematicReconstruction& b){
                             return  b.ntags < a.ntags || (b.ntags == a.ntags && b.weight < a.weight);
                         });
        
        sol_ = sols_[0];
    }
}



int KinematicReconstruction::getNSol()const
{
    return nSol_;
}



Struct_KinematicReconstruction KinematicReconstruction::getSol()const
{
    return sol_;
}



std::vector< Struct_KinematicReconstruction > KinematicReconstruction::getSols() const
{
    return sols_;
}



void KinematicReconstruction::loadData()
{
    std::cout<<"Smearing requires input distributions from files\n";
    
    r3_ = new TRandom3();
    
    // W mass
    TString data_path = common::DATA_PATH_COMMON();
    data_path.Append("/KinReco_wmass.root");
    TFile wmassfile(data_path);
    h_wmass_ = (TH1F*)wmassfile.Get("W_mass");
    h_wmass_->SetDirectory(0);
    wmassfile.Close();
    
    // jet and lepton resolution
    TString data_path1 = common::DATA_PATH_COMMON();
    data_path1.Append("/KinReco_d_angle_fE.root");
    TFile tempTFile(data_path1);
    h_jetAngleRes_ = (TH1F*)tempTFile.Get("KinReco_d_angle_lep_step7");
    h_jetEres_ = (TH1F*)tempTFile.Get("KinReco_fE_jet_step7");
    h_lepAngleRes_ = (TH1F*)tempTFile.Get("KinReco_d_angle_lep_step7");
    h_lepEres_ = (TH1F*)tempTFile.Get("KinReco_fE_lep_step7");
    h_lepAngleRes_->SetDirectory(0);
    h_lepEres_->SetDirectory(0);
    h_jetAngleRes_->SetDirectory(0);
    h_jetEres_->SetDirectory(0);
    tempTFile.Close();
     
// mbl
    TString data_path2 = common::DATA_PATH_COMMON();
        data_path2.Append("/KinReco_mbl.root");
    TFile fmbl(data_path2);
        h_mbl_w_ = (TH1F*)fmbl.Get("mbl_true");
        h_mbl_w_->SetDirectory(0);
        fmbl.Close();
        
    TString data_path21 = common::DATA_PATH_COMMON();
        data_path21.Append("/KinReco_mbl_wrong.root");
    TFile fmbl_wrong(data_path21);
        //h_mbl_w_ = (TH1F*)fmbl_wrong.Get("mbl_true_wrong");
        //h_mbl_w_->SetDirectory(0);
        fmbl_wrong.Close();
// ...
    
// pdf
        mttBinsForX_ = {340,400,470,550,650,750,850,950,1050,1150,1250,1550,1850,2150};
        TString data_path10 = common::DATA_PATH_COMMON();
        data_path10.Append("/KinReco_HERAPDF.root");
        TFile fpdf(data_path10);
        for(int i=0;i<(int)mttBinsForX_.size()-1;++i)
        {
            char name[20];
            sprintf(name,"h%d",(int)(mttBinsForX_.at(i+1)+mttBinsForX_.at(i))/2);
            h_pdf_w_[i] = (TH1F*)fpdf.Get(name);
            h_pdf_w_[i]->SetDirectory(0);
        }
        fpdf.Close();
// ...
    std::cout<<"Found all histograms needed for smearing\n";
}



void KinematicReconstruction::kinRecoMassLoop(const LV& leptonMinus, const LV& leptonPlus, const VLV* jets, const std::vector<double>* btags, const LV* met)
{

    std::vector<Struct_KinematicReconstruction> vect_sol;


    const TLorentzVector leptonPlus_tlv = common::LVtoTLV(leptonPlus);
    const TLorentzVector leptonMinus_tlv = common::LVtoTLV(leptonMinus);
    const TLorentzVector met_tlv = common::LVtoTLV(*met);

    std::vector<TLorentzVector> jets_tlv;
    for (const auto& jet : *jets) {
        jets_tlv.push_back(common::LVtoTLV(jet));
    }

    std::vector<int> b1_id;
    std::vector<int> b2_id;
    std::vector<double> btag_ww;
    std::vector<int> nb_tag;

    for(int i=0; i<(int)btags->size(); ++i) {
        for(int j=0; j<(int)btags->size(); ++j) {
            double wi = btags->at(i);
            double wj = btags->at(j);
//             if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0))continue;
            if(i==j || (wi<0.244 && wj<0.244))continue;
            btag_ww.push_back(wi + wj);

            if(wi>0.244 && wj>0.244){nb_tag.push_back(2); }
            else{nb_tag.push_back(1); }

            b1_id.push_back(i);
            b2_id.push_back(j);
        }
    }

    if(b1_id.size()<2) {
        nSol_=0;
        return;
    }

    for(int i=0; i<(int)btag_ww.size() - 1; ++i) {
        if(btag_ww[i]>=btag_ww[i+1]) continue;

        double aux = btag_ww[i];
        btag_ww[i] = btag_ww[i+1];
        btag_ww[i+1] = aux;
        int aix = b1_id[i];
        b1_id[i] = b1_id[i+1];
        b1_id[i+1] = aix;
        aix = b2_id[i];
        b2_id[i] = b2_id[i+1];
        b2_id[i+1] = aix;
        aix = nb_tag[i];
        nb_tag[i] = nb_tag[i+1];
        nb_tag[i+1] = aix;

        i=-1;

    }


    //jets loop

    nSol_=0;
    for(int ib=0; ib<(int)b1_id.size(); ++ib) {
        int j1=b1_id[ib];
        int j2=b2_id[ib];
//         if((*jets)[j1].Pt()<30 || (*jets)[j2].Pt()<30)continue;
//         if(fabs((*jets)[j1].Eta())>2.4 || fabs((*jets)[j2].Eta())>2.4)continue;
        const TLorentzVector l_temp=leptonMinus_tlv;
        const TLorentzVector al_temp=leptonPlus_tlv;
        const TLorentzVector b_temp=jets_tlv.at(j1);
        const TLorentzVector bbar_temp=jets_tlv.at(j2);
        const TLorentzVector met_temp=common::LVtoTLV(*met);
//         if((al_temp + b_temp).M()>180 || (l_temp + bbar_temp).M()>180)continue;



        // mass scan
        double vw_max = 0.;
        if(massLoop_){
           for(double iTopMass = 100.; iTopMass < 300.5; iTopMass += 1.){
//                 if(1) {double im=172.5;

                KinematicReconstruction_LSroutines tp_m(iTopMass, 4.8, 80.4, 0.0, 0.0);

                tp_m.setConstraints(al_temp, l_temp, b_temp, bbar_temp, met_temp.Px(), met_temp.Py());

//                 if(tp_m.GetNsol()<1 || tp_m.GetNsol()==1 || tp_m.GetNsol()==3) continue;
                if(tp_m.getNsol()<1) continue;
                
                if(!(tp_m.getTtSol()->at(0).weight>vw_max)) continue;
                
                nSol_++;
                
                vw_max=tp_m.getTtSol()->at(0).weight;
                sol_.jetB = b_temp;
                sol_.jetBbar = bbar_temp;
                sol_.lm = leptonMinus_tlv;
                sol_.lp = leptonPlus_tlv;
                sol_.met = met_temp;
                sol_.neutrino = tp_m.getTtSol()->at(0).neutrino;
                sol_.neutrinoBar = tp_m.getTtSol()->at(0).neutrinobar;
                sol_.weight = tp_m.getTtSol()->at(0).weight;
                sol_.Wplus = sol_.lp + sol_.neutrino;
                sol_.Wminus = sol_.lm + sol_.neutrinoBar;
                sol_.top = sol_.Wplus + sol_.jetB;
                sol_.topBar = sol_.Wminus + sol_.jetBbar;
                sol_.ttbar = sol_.top + sol_.topBar;
                sol_.jetB_index = j1;
                sol_.jetBbar_index = j2;
                sol_.ntags = nb_tag[ib];
            }
        }
        
        if(vw_max>0){
            vect_sol.push_back(sol_);
        }

    }

    if(nSol_>0){
        std::nth_element(begin(vect_sol), begin(vect_sol), end(vect_sol),
                         [](const Struct_KinematicReconstruction& a, const Struct_KinematicReconstruction& b){
                             return  b.ntags < a.ntags || (b.ntags == a.ntags && b.weight < a.weight);
                         });

        sol_=vect_sol[0];
    }

}






// -------------------------------------- Methods for KinematicReconstructionScaleFactors --------------------------------------

KinematicReconstructionScaleFactors::KinematicReconstructionScaleFactors(const std::vector<Channel::Channel>& channels,
                                                                         const Systematic::Systematic& systematic):
scaleFactor_(-999.)
{
    std::cout<<"--- Beginning preparation of kinematic reconstruction scale factors\n";
    
    // Set up proper internal systematic
    SystematicInternal systematicInternal(nominal);
    if(systematic.type() == Systematic::kin){
        if(systematic.variation() == Systematic::up) systematicInternal = vary_up;
        else if(systematic.variation() == Systematic::down) systematicInternal = vary_down;
    }
    
    // Set the scale factors according to specific systematic, and check whether all requested channels are defined
    this->prepareSF(systematicInternal);
    for(const auto& channel : channels){
        if(m_scaleFactor_.find(channel) == m_scaleFactor_.end()){
            std::cerr<<"ERROR in constructor of KinematicReconstructionScaleFactors! No scale factors defined for given channel: "
                     <<Channel::convert(channel)<<"\n...break\n"<<std::endl;
            exit(857);
        }
    }
    std::cout<<"Found scale factors for all requested channels\n";
    
    std::cout<<"=== Finishing preparation of kinematic reconstruction scale factors\n\n";
}



void KinematicReconstructionScaleFactors::prepareSF(const SystematicInternal& systematic)
{
    // FIXME: make proper documentation for SF determination (where?)
    // --> uncomment the following line to determine the Kin Reco SFs
    // --> then make && ./runNominalParallel.sh && ./Histo -t cp -p akr bkr step && ./kinRecoEfficienciesAndSF

    // SF = 1
    //const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 1.}, {Channel::emu, 1.}, {Channel::mumu, 1.} };
    //const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.}, {Channel::emu, 0.}, {Channel::mumu, 0.} };
    
    //SF for mass(top) = 100..300 GeV
    //const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 0.9779}, {Channel::emu, 0.9871}, {Channel::mumu, 0.9879} };
    //const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.0066}, {Channel::emu, 0.0032}, {Channel::mumu, 0.0056} };
    
    // SF for newKinReco flat Old
    //const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 0.9876}, {Channel::emu, 0.9921}, {Channel::mumu, 0.9949} };
    //const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.0043}, {Channel::emu, 0.0019}, {Channel::mumu, 0.0037} };
    
    // SF for newKinReco flat N007mvamet
    const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 0.9854}, {Channel::emu, 0.9934}, {Channel::mumu, 0.9934} };
    const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.0041}, {Channel::emu, 0.0018}, {Channel::mumu, 0.0036} };
    
    // SF for mass(top) = 173 GeV
    //const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 0.9696}, {Channel::emu, 0.9732}, {Channel::mumu, 0.9930} };
    //const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.0123}, {Channel::emu, 0.0060}, {Channel::mumu, 0.0105} };
    
    for(const auto& sfNominal : m_sfNominal){
        const Channel::Channel& channel = sfNominal.first;
        const double& centralValue = sfNominal.second;
        const double& uncertainty = m_sfUnc.at(channel);
        if(systematic == nominal) m_scaleFactor_[channel] = centralValue;
        else if(systematic == vary_up) m_scaleFactor_[channel] = centralValue + uncertainty;
        else if(systematic == vary_down) m_scaleFactor_[channel] = centralValue - uncertainty;
        else{
           std::cerr<<"ERROR in KinematicReconstructionScaleFactors::prepareSF()! "
                    <<"No treatment defined for specified systematic\n...break\n"<<std::endl;
           exit(857);
        }
    }
}



void KinematicReconstructionScaleFactors::prepareChannel(const Channel::Channel& channel)
{
    scaleFactor_ = m_scaleFactor_.at(channel);
}



double KinematicReconstructionScaleFactors::getSF()const
{
    return scaleFactor_;
}













