#include <cmath>
#include <cassert>
#include <algorithm>
#include <utility>
#include <cmath>

//#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>
// #include <TAttLine.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom3.h>

#include "classes.h"
#include "utils.h"
#include "analysisUtils.h"
// #include "KinReco.h"
#include "KinematicReconstruction.h"
#include "KinematicReconstruction_LSroutines.h"
#include "KinematicReconstruction_MeanSol.h"

#define printout 0

#define topmass (172.5)

using namespace std;









// -------------------------------------- Methods for KinematicReconstruction --------------------------------------



void KinematicReconstruction::angle_rot(double alpha, double e, TLorentzVector jet, TLorentzVector & jet_sm)
{
    /*Ganna Dolinska*/
  double px_1, py_1, pz_1; // sistema koordinat, gde impul's vdol' Oz
  
  /// Transition matrix detector -> syst1 ///
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  ///
  
  if(fabs(jet.Px())<=e){jet.SetPx(0);}
  if(fabs(jet.Py())<=e){jet.SetPy(0);}
  if(fabs(jet.Pz())<=e){jet.SetPz(0);}
  
  /// Rotation in syst 1 ///
//  int seed_random = (int)(( jet.Pt() - (int)(jet.Pt()) )*1000000000);
//  
//   TRandom3 r(seed_random);  ///random seed
  
  double phi = 2*TMath::Pi()*r3_->Rndm();
  pz_1 = jet.Vect().Mag()*cos(alpha);
  px_1 = - jet.Vect().Mag()*sin(alpha)*sin(phi);
  py_1 = jet.Vect().Mag()*sin(alpha)*cos(phi);  
  ///
  
  /// Transition detector <- syst1 ///
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
  ///
}


KinematicReconstruction::KinematicReconstruction():
nSol_(0),
isJetsMerging_(false),
h_wmass_(0),
h_jetAngleRes_(0),
h_jetEres_(0),
h_lepAngleRes_(0),
h_lepEres_(0),
//h_metAngleRes_(0),
//h_metPtres_(0),
//h_metPxRes_(0),
//h_metPyRes_(0),
h_nwcuts_(0),
//hvE_(0),
h_mbl_w_(0),
h_costheta_w_(0),
h_neuEta_w_(0)
{
    std::cout<<"--- Beginning preparation of kinematic reconstruction\n";
    this->loadData();
    std::cout<<"=== Finishing preparation of kinematic reconstruction\n\n";
}

void KinematicReconstruction::doJetsMerging(const VLV* jets,const std::vector<double> *btags)
{
    isJetsMerging_=true;   
    
    alljets_.clear();
    allbtags_.clear();
        for (const auto& jet : *jets) {
        alljets_.push_back(common::LVtoTLV(jet));
    }
    for(int i=0; i<(int)btags->size(); i++) {
        allbtags_.push_back(btags->at(i));
    }
    
}

void KinematicReconstruction::kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV *jets, const std::vector<double> *btags, const LV* met)
{
    
    sols_.clear();

    TLorentzVector leptonPlus_tlv = common::LVtoTLV(leptonPlus);
    TLorentzVector leptonMinus_tlv = common::LVtoTLV(leptonMinus);
    TLorentzVector met_tlv = common::LVtoTLV(*met);

        gRandom->SetSeed((int)(( leptonPlus_tlv.Pt() - (int)(leptonPlus_tlv.Pt()) )*1000000000));
    r3_->SetSeed(( leptonPlus_tlv.Pt() - (int)(leptonPlus_tlv.Pt()) )*1000000000);
    
    //jets selection

    vector<int> b1_id;
    vector<int> b2_id;
    vector<int> nb_tag;
    std::vector<TLorentzVector> jets_tlv;
    std::vector<double> new_btags;
    
    double btag_wp=0.244;
    
if(!isJetsMerging_)
{
            
    for (const auto& jet : *jets) {
        jets_tlv.push_back(common::LVtoTLV(jet));
    }
    
    for(int i=0; i<(int)btags->size(); i++) {
//         new_btags.push_back(btags->at(i));
        
        for(int j=0; j<(int)btags->size(); j++) {
            double wi = btags->at(i);
            double wj = btags->at(j);
//                  if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0))continue;
            if(i==j || (wi<btag_wp && wj<btag_wp))continue;

            if(wi>btag_wp && wj>btag_wp){nb_tag.push_back(2); }
            else{nb_tag.push_back(1); }

            b1_id.push_back(i);
            b2_id.push_back(j);
        }
    }
}
else
{
    
    //fill new jets:  logik I.
//     std::vector<TLorentzVector> new_jets_tlv;
//     std::vector<double> new_btag;
//     std::vector<int> veto_i;
//     std::vector<int> veto_j;
//         for(int i=0; i<(int)alljets_.size(); i++)
//         {
//             if(allbtags_.at(i)<btag_wp)continue;
//             TLorentzVector ijet=alljets_.at(i);
//             
//             for(int j=0; j<(int)alljets_.size(); j++)
//             {
//                  if(i==j)continue;
//                  if(allbtags_.at(j)>=btag_wp)continue; 
//                 TLorentzVector jjet=alljets_.at(j);
//                 double dRij=ijet.DeltaR(jjet);
//                  if(dRij>1)continue;
//                 TLorentzVector new_jet=ijet+jjet;
//                     if(new_jet.Pt()<30)continue;
//                     if(fabs(new_jet.Eta())>2.4)continue;
//                 new_jets_tlv.push_back(new_jet);
//                 new_btag.push_back(allbtags_.at(i));
//                 veto_i.push_back(i);
//                 veto_j.push_back(j);
//             }   
//         }
//  
//      std::vector<int> index_i;
//      std::vector<int> index_j;
//      for(int i=0; i<(int)alljets_.size(); i++)
//         {
//                 if(alljets_.at(i).Pt()<30)continue;
//                 if(fabs(alljets_.at(i).Eta())>2.4)continue;
//                 jets_tlv.push_back(alljets_.at(i));
//                 new_btags.push_back(allbtags_.at(i));
//                 index_i.push_back(i);
//                 index_j.push_back(i);
//         }
//     
//     
//     for(int i=0; i<(int)new_jets_tlv.size(); i++)
//     {
//         jets_tlv.push_back(new_jets_tlv.at(i));
//         new_btags.push_back(new_btag.at(i));
//         index_i.push_back(veto_i.at(i));
//         index_j.push_back(veto_j.at(i));
//     }
//     
//         for(int i=0; i<(int)jets_tlv.size(); i++) 
//         {
//                 for(int j=0; j<(int)jets_tlv.size(); j++) 
//                 {
//                                 if(index_i[i]==index_i[j] || index_i[i]==index_j[j] || index_j[i]==index_i[j] || index_j[i]==index_j[j])continue;                   
//                                double wi = new_btags.at(i);
//                                double wj = new_btags.at(j);
// //                                 //if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0))continue;
//                                 if(i==j || (wi<btag_wp && wj<btag_wp))continue;
//                                 if(wi>btag_wp && wj>btag_wp){nb_tag.push_back(2); }
//                                 else{nb_tag.push_back(1); }
//                                 b1_id.push_back(i);
//                                 b2_id.push_back(j);
//                 }
//         }


//fill new jets:  logik II.

    std::vector<TLorentzVector> new_jets_tlv;
    std::vector<double> new_btag;
    std::vector<int> index_dR;
    int index_i;
    int index_j;
    
    for(int i=0; i<(int)alljets_.size(); i++)
    {
            index_dR.push_back(-1);

    }
    
   double min_dR=1;
   do{
       min_dR=1;
       index_i=-1;
       index_j=-1;
        for(int i=0; i<(int)alljets_.size(); i++)
        {
            if(alljets_.at(i).Pt()<10)continue;
            if(index_dR[i]>-1)continue;
            TLorentzVector ijet=alljets_.at(i);
             
            for(int j=0; j<(int)alljets_.size(); j++)
            {
                if(i==j)continue;
                if(alljets_.at(j).Pt()<10)continue;
                if(index_dR[j]>-1)continue;
                TLorentzVector jjet=alljets_.at(j);
                
                double dRij=ijet.DeltaR(jjet);
                
                if(dRij>=1)continue;
                
                if(dRij<min_dR)
                {
                    min_dR=dRij;
                    index_i=i;
                    index_j=j;
                    
                } 
            }
        }
        if(index_i>=0)
        {
            index_dR[index_i]=index_j;
            index_dR[index_j]=index_i;
            index_i=-1;
            index_j=-1;
        }
        
      }   
     while(min_dR<1);
    
     
     vector<int> used_index;
     
     for(int i=0; i<(int)alljets_.size(); i++)
     {
         bool flag=false;
         for(int k=0; k<(int)used_index.size(); k++)
         {
             if(used_index[k]==i)flag=true;
        }
        if(flag)continue;
         
        if(index_dR[i]<0)
        {
                if(alljets_.at(i).Pt()<30)continue;
                if(fabs(alljets_.at(i).Eta())>2.4)continue;
                jets_tlv.push_back(alljets_.at(i));
                new_btags.push_back(allbtags_.at(i));
                
        }
        else
        {
            if(index_dR[i]>-1)
            {
                TLorentzVector temp = alljets_.at(i) + alljets_.at((index_dR[i]));
                used_index.push_back(index_dR[i]);
                if(temp.Pt()<30)continue;
                if(fabs(temp.Eta())>2.4)continue;
                jets_tlv.push_back(temp);
                new_btags.push_back(allbtags_.at(i)*(allbtags_.at(i)>allbtags_.at(index_dR[i])) + allbtags_.at(index_dR[i])*(allbtags_.at(index_dR[i])>allbtags_.at(i)));
            }
        }
        
        
     }
     
        for(int i=0; i<(int)jets_tlv.size(); i++) 
        {
                for(int j=0; j<(int)jets_tlv.size(); j++) 
                {
                               double wi = new_btags.at(i);
                               double wj = new_btags.at(j);
//                                 //if(i==j || (wi<0.244 && wj<0.244) || (wi<0 || wj<0))continue;
                                if(i==j || (wi<btag_wp && wj<btag_wp))continue;
                                if(wi>btag_wp && wj>btag_wp){nb_tag.push_back(2); }
                                else{nb_tag.push_back(1); }
                                b1_id.push_back(i);
                                b2_id.push_back(j);
                }
        }
     
    
} //else !isJetsMerging_


if(b1_id.size()<2)return;     

    KinematicReconstruction_MeanSol meanSol(topmass);
    //double max_sum_weight=0;
    int /*nbtag=0,*/ isHaveSol=0;

    for(int ib=0; ib<(int)b1_id.size(); ib++) {
        isHaveSol=0;
        int j1=b1_id[ib];
        int j2=b2_id[ib];
//         if((*jets)[j1].Pt()<30 || (*jets)[j2].Pt()<30)continue;
//         if(fabs((*jets)[j1].Eta())>2.4 || fabs((*jets)[j2].Eta())>2.4)continue;
        TLorentzVector l_temp, al_temp, b_temp, bbar_temp, met_temp;
        l_temp=leptonMinus_tlv;
        al_temp=leptonPlus_tlv;
        b_temp=jets_tlv.at(j1);
        bbar_temp=jets_tlv.at(j2);
        met_temp.SetXYZM(met->Px(), met->Py(), 0, 0);

        if((al_temp + b_temp).M()>180 || (l_temp + bbar_temp).M()>180)continue;

//           double b_w=new_btags.at(j1)/(1-new_btags.at(j1));
//           double bbar_w=new_btags.at(j2)/(1-new_btags.at(j2));
//             double b_w=1;
//          double bbar_w=1;


        TVector3 vX_reco =  - b_temp.Vect() - bbar_temp.Vect() - l_temp.Vect() - al_temp.Vect() - met_temp.Vect();       

        for(int sm=0; sm<100; sm++) {
            TLorentzVector b_sm=b_temp;
            TLorentzVector bbar_sm=bbar_temp;
            TLorentzVector met_sm;//=met_temp;
            TLorentzVector l_sm=l_temp;
            TLorentzVector al_sm=al_temp;
            TLorentzVector vX_sm;

            double fB=h_jetEres_->GetRandom();//fB=1;  //sm off
            double xB=sqrt((fB*fB*b_sm.E()*b_sm.E()-b_sm.M2())/(b_sm.P()*b_sm.P()));
            double fBbar=h_jetEres_->GetRandom();//fBbar=1; //sm off
            double xBbar=sqrt((fBbar*fBbar*bbar_sm.E()*bbar_sm.E()-bbar_sm.M2())/(bbar_sm.P()*bbar_sm.P()));
                           
            double fL=h_lepEres_->GetRandom();//fL=1; //sm off
            double xL=sqrt((fL*fL*l_sm.E()*l_sm.E()-l_sm.M2())/(l_sm.P()*l_sm.P()));
            double faL=h_lepEres_->GetRandom();//faL=1;  //sm off
            double xaL=sqrt((faL*faL*al_sm.E()*al_sm.E()-al_sm.M2())/(al_sm.P()*al_sm.P()));
                           
            b_sm.SetXYZT(b_sm.Px()*xB,b_sm.Py()*xB,b_sm.Pz()*xB,b_sm.E()*fB);
            angle_rot(h_jetAngleRes_->GetRandom(),0.001,b_sm,b_sm);
                           
            bbar_sm.SetXYZT(bbar_sm.Px()*xBbar,bbar_sm.Py()*xBbar,bbar_sm.Pz()*xBbar,bbar_sm.E()*fBbar);    
            angle_rot(h_jetAngleRes_->GetRandom(),0.001,bbar_sm,bbar_sm);
                                
            l_sm.SetXYZT(l_sm.Px()*xL,l_sm.Py()*xL,l_sm.Pz()*xL,l_sm.E()*fL);
            angle_rot(h_lepAngleRes_->GetRandom(),0.001,l_sm,l_sm);
                           
            al_sm.SetXYZT(al_sm.Px()*xaL,al_sm.Py()*xaL,al_sm.Pz()*xaL,al_sm.E()*faL);
            angle_rot(h_lepAngleRes_->GetRandom(),0.001,al_sm,al_sm);
                           
// met Pt+angle smearing
//             double fX=1;
//             double dAX=0;
//             for(int i=0;i<((int)(ptBins_.size()));i++)
//             {
//                   if(vX_reco.Pt()>ptBins_[i]&&vX_reco.Pt()<=ptBins_[i+1])
//                   {
//                        //fX=h_metPtres_[i]->GetRandom();
//                        //dAX=h_metAngleRes_[i]->GetRandom();
//                        
//                        fX=1;
//                         dAX=0;
//                   }
//             }
//             double tempX= vX_reco.Px()*fX;
//             double tempY= vX_reco.Py()*fX;
//             if(r.Rndm()>0.5)dAX=-dAX;
//             vX_sm.SetXYZM(tempX*cos(dAX)-tempY*sin(dAX),tempX*sin(dAX)+tempY*cos(dAX),0,0);
//                            
// ...

//met Pt+angle smearing
            double fPx=1;
            double fPy=0;
            for(int i=0;i<((int)(ptBins_.size()));i++)
            {
                  if(vX_reco.Pt()>ptBins_[i]&&vX_reco.Pt()<=ptBins_[i+1])
                  {
                       //fPx=h_metPxRes_[i]->GetRandom();
                       //fPy=h_metPyRes_[i]->GetRandom();
                       fPx=1;
                       fPy=1;
                  }
            }
            double tempX= vX_reco.Px()*fPx;
            double tempY= vX_reco.Py()*fPy;
            vX_sm.SetXYZM(tempX,tempY,0,0);
       
//...
            TVector3 metV3_sm= -b_sm.Vect()-bbar_sm.Vect()-l_sm.Vect()-al_sm.Vect()-vX_sm.Vect();
            met_sm.SetXYZM(metV3_sm.Px(),metV3_sm.Py(),0,0);
            
            
            KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom());
//            KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom(),(*h_nwcuts_));
            //KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom(),(*h_neuEta_w_));
//              KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom(),hvE_);
//             KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom(),hvE_,(*h_neuEta_w_));
//              KinematicReconstruction_LSroutines tp_sm(0.0,0.0,h_wmass_->GetRandom(),h_wmass_->GetRandom(),(*h_costheta_w_),1);
            
            tp_sm.setConstraints(al_sm, l_sm, b_sm, bbar_sm, met_sm.Px(), met_sm.Py());


            if(tp_sm.getNsol()>0)
            {
                isHaveSol=1;
                for(int i=0; i<=tp_sm.getNsol()*0; i++) 
                {
                   // meanSol.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,1);
                  //meanSol.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,(tp_sm.GetTtSol()->at(i).lepEw));//*b_w*bbar_w
                  double mbl_weight=h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((al_sm+b_sm).M()))*h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((l_sm+bbar_sm).M()))/100000000;
                  meanSol.add(tp_sm.getTtSol()->at(i).top,tp_sm.getTtSol()->at(i).topbar,tp_sm.getTtSol()->at(i).neutrino,tp_sm.getTtSol()->at(i).neutrinobar,mbl_weight);//
                  //meanSol.add(tp_sm.GetTtSol()->at(i).top,tp_sm.GetTtSol()->at(i).topbar,tp_sm.GetTtSol()->at(i).neutrino,tp_sm.GetTtSol()->at(i).neutrinobar,(tp_sm.GetTtSol()->at(i).vw));//

                   
                }
            }
            
			} // for sm=100 end

//         if(isHaveSol) {
//             
//                 if(nb_tag[ib]==nbtag){
//                     if(meanSol.getSumWeight()>max_sum_weight) {
//                         max_sum_weight=meanSol.getSumWeight();
//                         meanSol.getMeanSol(sol_.top, sol_.topBar, sol_.neutrino, sol_.neutrinoBar);
//                         sol_.jetB = b_temp;
//                         sol_.jetBbar = bbar_temp;
//                         sol_.lm = leptonMinus_tlv;
//                         sol_.lp = leptonPlus_tlv;
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
//                 sol_.lm = leptonMinus_tlv;
//                 sol_.lp = leptonPlus_tlv;
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

        if(isHaveSol) 
            {
                        meanSol.getMeanSol(sol_.top, sol_.topBar, sol_.neutrino, sol_.neutrinoBar);
                        sol_.jetB = b_temp;
                        sol_.jetBbar = bbar_temp;
                        sol_.lm = leptonMinus_tlv;
                        sol_.lp = leptonPlus_tlv;
                        sol_.met = met_temp;
                        sol_.weight = meanSol.getSumWeight();
                        sol_.Wplus = sol_.lp + sol_.neutrino;
                        sol_.Wminus = sol_.lm + sol_.neutrinoBar;
                        sol_.ttbar = sol_.top + sol_.topBar;
                        sol_.jetB_index = b1_id[ib];
                        sol_.jetBbar_index = b2_id[ib];
                        sol_.ntags = nb_tag[ib];
                sols_.push_back(sol_);
            }


        meanSol.clear();

    } ///end jets loop

    nSol_=(int)(sols_.size());
    
     if(nSol_>0){
        std::nth_element(begin(sols_), begin(sols_), end(sols_),
                         [](const Struct_KinematicReconstruction& a, const Struct_KinematicReconstruction& b){
                             return  b.ntags < a.ntags || (b.ntags == a.ntags && b.weight < a.weight);
                         });

        sol_=sols_[0];
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

vector< Struct_KinematicReconstruction > KinematicReconstruction::getSols() const
{
    return sols_;
}


void KinematicReconstruction::loadData()
{
    
    r3_ = new TRandom3();
    
//     // W mass

    TString data_path = common::DATA_PATH_COMMON();
    data_path.Append("/KinReco_wmass.root");
    TFile wmassfile(data_path);
    h_wmass_ = (TH1F*)wmassfile.Get("W_mass");
    h_wmass_->SetDirectory(0);
    wmassfile.Close();
    
    // jet resolution
    TString data_path1 = common::DATA_PATH_COMMON();
    data_path1.Append("/KinReco_d_angle_jet.root");
    TFile danglejetfile(data_path1);
    h_jetAngleRes_ = (TH1F*)danglejetfile.Get("d_angle_jet");
    h_jetEres_ = (TH1F*)danglejetfile.Get("fE");
    h_jetAngleRes_->SetDirectory(0);
    h_jetEres_->SetDirectory(0);
    danglejetfile.Close();
    //
    
    //lepton resolution
    TString data_path2 = common::DATA_PATH_COMMON();
    data_path2.Append("/KinReco_d_angle_lep.root");
    TFile lepangleresfile(data_path2);
    h_lepAngleRes_ = (TH1F*)lepangleresfile.Get("d_angle_jet");
    h_lepEres_ = (TH1F*)lepangleresfile.Get("fE");
    h_lepAngleRes_->SetDirectory(0);
    h_lepEres_->SetDirectory(0);
    lepangleresfile.Close();
    //
// 
//     //MET resolution
    
    ptBins_ =  {0,7,13,18,27,37,48,59,70,100,200,500,1000,1500};
    TString data_path3 = common::DATA_PATH_COMMON();
    data_path3.Append("/KinReco_rmshistA_dAngle_vs_Pt.root");
    TFile metangleresfile(data_path3);
    for(int i=0;i<13;i++)
    {
        char temp_name[100];
        sprintf(temp_name,"bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metAngleRes_[i]=(TH1F*)metangleresfile.Get(temp_name);
        sprintf(temp_name,"metAng_bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metAngleRes_[i]->SetName(temp_name);
        h_metAngleRes_[i]->SetDirectory(0);
    }
    metangleresfile.Close();
    
    TString data_path4 = common::DATA_PATH_COMMON();
    data_path4.Append("/KinReco_rmshistA_fPt_vs_Pt.root");
    TFile metptresfile(data_path4);
    for(int i=0;i<13;i++)
    {
        char temp_name[100];
        sprintf(temp_name,"bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPtres_[i]=(TH1F*)metptresfile.Get(temp_name);
        sprintf(temp_name,"metPt_bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPtres_[i]->SetName(temp_name);
        h_metPtres_[i]->SetDirectory(0);
    }
    metptresfile.Close();
    
   //met px res
    TString data_path41 = common::DATA_PATH_COMMON(); 
    data_path41.Append("/KinReco_rmshist_fPx_vs_Pt.root");
    TFile metpxresfile(data_path41);
    for(int i=0;i<13;i++)
    {
        char temp_name[100];
        sprintf(temp_name,"bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPxRes_[i]=(TH1F*)metpxresfile.Get(temp_name);
        sprintf(temp_name,"metPx_bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPxRes_[i]->SetName(temp_name);
        h_metPxRes_[i]->SetDirectory(0);
    }
    metpxresfile.Close();
  //...   
    
  //met py res
    TString data_path42 = common::DATA_PATH_COMMON(); 
    data_path42.Append("/KinReco_rmshist_fPy_vs_Pt.root");
    TFile metpyresfile(data_path42);
    for(int i=0;i<13;i++)
    {
        char temp_name[100];
        sprintf(temp_name,"bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPyRes_[i]=(TH1F*)metpyresfile.Get(temp_name);
        sprintf(temp_name,"metPy_bin_%.0f_%.0f",ptBins_[i],ptBins_[i+1]);
        h_metPyRes_[i]->SetName(temp_name);
        h_metPyRes_[i]->SetDirectory(0);
    }
    metpyresfile.Close();
  //...   


    TString data_path5 = common::DATA_PATH_COMMON();
    data_path5.Append("/KinReco_wneutrinocuts.root");
    TFile fnw(data_path5);
    h_nwcuts_ = (TH1F*)fnw.Get("Eneu_true");
    h_nwcuts_->SetDirectory(0);
    fnw.Close();
    
///  E 1d bins
        TString data_path6 = common::DATA_PATH_COMMON();
        data_path6.Append("/KinReco_wneutrino_E.root"); // no norm
//         data_path6.Append("/KinReco_wneutrino_E_normilyze.root"); // norm
        TFile frootNeut(data_path6);
        char hvE_name[20];
        
                        for(int i=0;i<6;i++)
                        {
                            sprintf(hvE_name,"Eneu_true%d",i+1);
                                    hvE_[i] = (TH1F*)frootNeut.Get(hvE_name);
                                    hvE_[i]->SetDirectory(0);
                        }
    frootNeut.Close();
///...

/// mbl

        TString data_path7 = common::DATA_PATH_COMMON();
        data_path7.Append("/KinReco_mbl.root");
        TFile fmbl(data_path7);
        h_mbl_w_ = (TH1F*)fmbl.Get("mbl_true");
        h_mbl_w_->SetDirectory(0);
        fmbl.Close();
        
        TString data_path71 = common::DATA_PATH_COMMON();
        data_path71.Append("/KinReco_mbl_wrong.root");
        TFile fmbl_wrong(data_path71);
//         h_mbl_w_ = (TH1F*)fmbl_wrong.Get("mbl_true_wrong");
//         h_mbl_w_->SetDirectory(0);
        fmbl_wrong.Close();
        
/// ...        
        
/// cos_theta*

        TString data_path8 = common::DATA_PATH_COMMON();
        data_path8.Append("/KinReco_costheta_true.root");
        TFile fcostheta(data_path8);
        h_costheta_w_ = (TH1F*)fcostheta.Get("cos_theta_true");
        h_costheta_w_->SetDirectory(0);
        fcostheta.Close();
/// ...
        
/// neuEta
        TString data_path9 = common::DATA_PATH_COMMON();
        data_path9.Append("/KinReco_wneutrinoEta.root");
        TFile fneuEta(data_path9);
        h_neuEta_w_ = (TH1F*)fneuEta.Get("Etaneu_true");
        h_neuEta_w_->SetDirectory(0);
        fneuEta.Close();
/// ...        
        
        
}

void KinematicReconstruction::kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV *jets, const std::vector<double> *btags, const LV* met, bool mass_loop_on)
{

    std::vector<Struct_KinematicReconstruction> vect_sol;


    TLorentzVector leptonPlus_tlv = common::LVtoTLV(leptonPlus);
    TLorentzVector leptonMinus_tlv = common::LVtoTLV(leptonMinus);
    TLorentzVector met_tlv = common::LVtoTLV(*met);

    std::vector<TLorentzVector> jets_tlv;
    for (const auto& jet : *jets) {
        jets_tlv.push_back(common::LVtoTLV(jet));
    }

    vector<int> b1_id;
    vector<int> b2_id;
    vector<double> btag_ww;
    vector<int> nb_tag;

    for(int i=0; i<(int)btags->size(); i++) {
        for(int j=0; j<(int)btags->size(); j++) {
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

    for(int i=0; i<(int)btag_ww.size() - 1; i++) {
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


    ///jets loop

    nSol_=0;
    for(int ib=0; ib<(int)b1_id.size(); ib++) {
        int j1=b1_id[ib];
        int j2=b2_id[ib];
//         if((*jets)[j1].Pt()<30 || (*jets)[j2].Pt()<30)continue;
//         if(fabs((*jets)[j1].Eta())>2.4 || fabs((*jets)[j2].Eta())>2.4)continue;
        TLorentzVector l_temp, al_temp, b_temp, bbar_temp, met_temp;
        l_temp=leptonMinus_tlv;
        al_temp=leptonPlus_tlv;
        b_temp=jets_tlv.at(j1);
        bbar_temp=jets_tlv.at(j2);
        met_temp.SetXYZM(met->Px(), met->Py(), 0, 0);
//         if((al_temp + b_temp).M()>180 || (l_temp + bbar_temp).M()>180)continue;

        /*smearing*/
        double vw_max=0;
        if(!mass_loop_on){

//             TF1 fJet("fJet", "[0]*pow((x - 0), 0.5) + [1]");
//             fJet.SetParameters(4.18432, - 24.8944);
//
//             TF1 fLep("fLep", "[0]*pow((x - 0), 0.5) + [1]");
//             fLep.SetParameters(0.575716, - 3.11214);


            TRandom3 r(0); ///random seed

//             TF1 f1("f1", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x", 30, 580);
//             f1.SetParameters(1.06517e + 01, 3.58131e - 02, 7.44249e - 04, - 1.70968e - 06, 1.21316e - 09);

            double bRMS=(4.13*pow((b_temp.E() - 0), 0.5) - 27.32)/b_temp.E(); //fJet.Eval(b_temp.E())/b_temp.E(); //%
            double bbarRMS=(4.13*pow((bbar_temp.E() - 0), 0.5) - 27.32)/bbar_temp.E(); //fJet.Eval(bbar_temp.E())/bbar_temp.E(); //%
            double lRMS=(0.575716*pow((l_temp.E() - 0), 0.5) - 3.11214)/l_temp.E(); //fLep.Eval(l_temp.E())/l_temp.E(); //%
            double alRMS=(0.575716*pow((al_temp.E() - 0), 0.5) - 3.11214)/al_temp.E(); //fLep.Eval(al_temp.E())/al_temp.E(); //%


            TVector3 vX_reco =  - b_temp.Vect() - bbar_temp.Vect() - l_temp.Vect() - al_temp.Vect() - met_temp.Vect();

            double xRMS=20/vX_reco.Pt(); //%
            for(int sm=0; sm<100; sm++) {
                TLorentzVector b_sm=b_temp;
                TLorentzVector bbar_sm=bbar_temp;
                TLorentzVector met_sm=met_temp;
                TLorentzVector l_sm=l_temp;
                TLorentzVector al_sm=al_temp;


                double fB=r.Gaus(1, bRMS);
                double xB=sqrt((fB*fB*b_sm.E()*b_sm.E() - b_sm.M2())/(b_sm.P()*b_sm.P()));
                double fBbar=r.Gaus(1, bbarRMS);
                double xBbar=sqrt((fBbar*fBbar*bbar_sm.E()*bbar_sm.E() - bbar_sm.M2())/(bbar_sm.P()*bbar_sm.P()));

                double fL=r.Gaus(1, lRMS);
                double xL=sqrt((fL*fL*l_sm.E()*l_sm.E() - l_sm.M2())/(l_sm.P()*l_sm.P()));
                double faL=r.Gaus(1, alRMS);
                double xaL=sqrt((faL*faL*al_sm.E()*al_sm.E() - al_sm.M2())/(al_sm.P()*al_sm.P()));

                double xX=r.Gaus(1, xRMS);

//                 met_sm.SetXYZM(((b_sm.Px()*(1 - xB) + bbar_sm.Px()*(1 - xBbar))/(met_sm.Px()) + 1)*met_sm.Px(), ((b_sm.Py()*(1 - xB) + bbar_sm.Py()*(1 - xBbar))/(met_sm.Py()) + 1)*met_sm.Py(), 0, 0); // tollko b - jet
//                 met_sm.SetXYZM(((b_sm.Px()*(1 - xB) + bbar_sm.Px()*(1 - xBbar) + l_sm.Px()*(1 - xL) + al_sm.Px()*(1 - xaL))/(met_sm.Px()) + 1)*met_sm.Px(), ((b_sm.Py()*(1 - xB) + bbar_sm.Py()*(1 - xBbar) + l_sm.Py()*(1 - xL) + al_sm.Py()*(1 - xaL))/(met_sm.Py()) + 1)*met_sm.Py(), 0, 0); //only b - jet and lep

                met_sm.SetXYZM(((b_sm.Px()*(1 - xB) + bbar_sm.Px()*(1 - xBbar) + l_sm.Px()*(1 - xL) + al_sm.Px()*(1 - xaL) + vX_reco.Px()*(1 - xX))/(met_sm.Px()) + 1)*met_sm.Px(), ((b_sm.Py()*(1 - xB) + bbar_sm.Py()*(1 - xBbar) + l_sm.Py()*(1 - xL) + al_sm.Py()*(1 - xaL) + vX_reco.Py()*(1 - xX))/(met_sm.Py()) + 1)*met_sm.Py(), 0, 0); //all

                b_sm.SetXYZT(b_sm.Px()*xB, b_sm.Py()*xB, b_sm.Pz()*xB, b_sm.E()*fB);
                bbar_sm.SetXYZT(bbar_sm.Px()*xBbar, bbar_sm.Py()*xBbar, bbar_sm.Pz()*xBbar, bbar_sm.E()*fBbar);

                l_sm.SetXYZT(l_sm.Px()*xL, l_sm.Py()*xL, l_sm.Pz()*xL, l_sm.E()*fL);
                al_sm.SetXYZT(al_sm.Px()*xaL, al_sm.Py()*xaL, al_sm.Pz()*xaL, al_sm.E()*faL);

                KinematicReconstruction_LSroutines tp_sm(topmass, 4.8, 80.4, 0.0, 0.0);
                tp_sm.setConstraints(al_sm, l_sm, b_sm, bbar_sm, met_sm.Px(), met_sm.Py());

                if(!(tp_sm.getNsol()<1 || tp_sm.getNsol()==1 || tp_sm.getNsol()==3)) {
                    if(tp_sm.getTtSol()->at(0).weight>vw_max) {
                        nSol_++;

                        vw_max=tp_sm.getTtSol()->at(0).weight;
                        sol_.jetB = b_temp;
                        sol_.jetBbar = bbar_temp;
                        sol_.lm = leptonMinus_tlv;
                        sol_.lp = leptonPlus_tlv;
                        sol_.met = met_temp;
                        sol_.neutrino = tp_sm.getTtSol()->at(0).neutrino;
                        sol_.neutrinoBar = tp_sm.getTtSol()->at(0).neutrinobar;
                        sol_.weight = tp_sm.getTtSol()->at(0).weight;
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
            }
        }


        /*smearing*/
        /*mass scan*/

        if(mass_loop_on){
           for(double im=100; im<300.5; im+=1){
//                 if(1) {double im=172.5;

                KinematicReconstruction_LSroutines tp_m(im, 4.8, 80.4, 0.0, 0.0);

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

        /*mass scan*/

        if(vw_max>0){
            vect_sol.push_back(sol_);
        }

    } ///end jets loop

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
    
    // SF for newKinReco flat
    const std::map<Channel::Channel, double> m_sfNominal { {Channel::ee, 0.9876}, {Channel::emu, 0.9921}, {Channel::mumu, 0.9949} };
    const std::map<Channel::Channel, double> m_sfUnc { {Channel::ee, 0.0043}, {Channel::emu, 0.0019}, {Channel::mumu, 0.0037} };
    
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













