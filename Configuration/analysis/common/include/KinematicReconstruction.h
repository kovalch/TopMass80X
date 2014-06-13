#ifndef KinematicReconstruction_h
#define KinematicReconstruction_h

#include <vector>
#include <map>
#include <iostream>
#include <stdio.h>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TH1F.h>

#include <Math/PtEtaPhiM4D.h>
#include <Math/LorentzVector.h>
#include <Rtypes.h>
#include <TRandom3.h>
#include "classesFwd.h"
#include "sampleHelpers.h"





struct Struct_KinematicReconstruction {
    TLorentzVector lp, lm;
    TLorentzVector jetB, jetBbar;
    size_t jetB_index, jetBbar_index;
    TLorentzVector met;
    TLorentzVector Wplus, Wminus;
    TLorentzVector top, topBar;
    TLorentzVector neutrino, neutrinoBar;
    TLorentzVector ttbar;
    double recMtop;
    double weight;
    int ntags;
};



class KinematicReconstruction {

public:

  KinematicReconstruction();
  
int getNSol()const;
Struct_KinematicReconstruction getSol()const;
std::vector<Struct_KinematicReconstruction> getSols()const;

void loadData();
void kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV *jets, const std::vector<double> *btags, const LV* met);
void kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV *jets, const std::vector<double> *btags, const LV* met, bool mass_loop_on);
void doJetsMerging(const VLV *jets,const std::vector<double> *btags); 

private:

    TRandom3* r3_;
    void angle_rot(double alpha, double e, TLorentzVector jet, TLorentzVector & jet_sm);
    
    int nSol_;
    Struct_KinematicReconstruction sol_;
    std::vector<Struct_KinematicReconstruction> sols_;
    
    bool isJetsMerging_;
    std::vector<TLorentzVector> alljets_;
    std::vector<double> allbtags_;
//     // W mass
    TH1F * h_wmass_;
// 
//     // jet resolution
    TH1F * h_jetAngleRes_;
    TH1F * h_jetEres_;
//        
//     //lepton resolution
    TH1F * h_lepAngleRes_;
    TH1F * h_lepEres_;
// 
//     //MET resolution
//     double ptBins[14];
    std::vector<double> ptBins_;
    TH1F * h_metAngleRes_[13];
    TH1F * h_metPtres_[13];
    
    TH1F * h_metPxRes_[13];
    TH1F * h_metPyRes_[13];
    
// 
    TH1F *h_nwcuts_;
// 
// //  E 1d bins
        TH1F* hvE_[6];
// 
// // mbl
        TH1F *h_mbl_w_;
// // costheta        
        TH1F *h_costheta_w_;

// neuEta 0d weight
        TH1F *h_neuEta_w_;
//
};







class KinematicReconstructionScaleFactors{
    
public:
    
    /// Constructor
    KinematicReconstructionScaleFactors(const std::vector<Channel::Channel>& channels,
                                        const Systematic::Systematic& systematic);
    
    /// Destructor
    ~KinematicReconstructionScaleFactors(){}
    
    
    
    /// Prepare scale factors per channel
    void prepareChannel(const Channel::Channel& channel);
    
    /// Get kinematic reconstruction per-event scale factor
    double getSF()const;
    
    
    
private:
    
    /// Enumeration for possible systematics
    enum SystematicInternal{nominal, vary_up, vary_down, undefined};
    
    /** Set up the scale factor for the Kinematic Reconstruction
     *
     * Currently a flat per-channel SF is used. For the systematic KIN_UP and KIN_DOWN,
     * the SF is modified by its uncertainty.
     *
     * To calculate the SF, you need to set the SF to 1 and rerun. Then determine the
     * SF with kinRecoEfficienciesAndSF
     */
    void prepareSF(const SystematicInternal& systematic);
    
    
    
    /// Map containing the flat scale factors for all channels
    std::map<Channel::Channel, double> m_scaleFactor_;
    
    /// The per-channel scale factor
    double scaleFactor_;
};








#endif
