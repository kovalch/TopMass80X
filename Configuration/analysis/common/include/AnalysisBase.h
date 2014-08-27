#ifndef AnalysisBase_h
#define AnalysisBase_h

#include <string>
#include <vector>
#include <functional>

#include <TSelector.h>
#include <Rtypes.h>
#include <TString.h>

class TBranch;
class TTree;
class TH1;

#include "classesFwd.h"
#include "ScaleFactorsFwd.h"
#include "storeTemplate.h"
#include "sampleHelpers.h"

class KinematicReconstruction;
class KinematicReconstructionScaleFactors;
class PileupScaleFactors;
class TriggerScaleFactors;
class LeptonScaleFactors;
class BtagScaleFactors;
class JetEnergyResolutionScaleFactors;
class JetEnergyScaleScaleFactors;
class TopPtScaleFactors;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class ZGenObjects;
class KinRecoObjects;
class MetRecoilCorrector;





class AnalysisBase : public TSelector
{
    
public:

    /// Constructor
    AnalysisBase(TTree* =0);

    /// Destructor
    virtual ~AnalysisBase(){};

    /// Inherited TSelector methods
    virtual void Begin(TTree*);
    virtual void SlaveBegin(TTree*);
    virtual Bool_t Process(Long64_t);
    virtual void Init(TTree*);
    virtual void SlaveTerminate();
    virtual void Terminate();

    /// Set sample name
    void SetSamplename(const TString& samplename);
    
    /// Set generator parameters
    void SetGeneratorBools(const TString& samplename, const Systematic::Systematic& systematic);
    
    /// Set Channel
    void SetChannel(const Channel::Channel& channel);
    
    /// Set systematic
    void SetSystematic(const Systematic::Systematic& systematic);
    
    /// Set whether it is MC sample
    void SetMC(const bool isMC);
    
    /// Set whether it is Top signal sample
    void SetTopSignal(const bool isTopSignal);
    
    /// Set whether it is Higgs signal sample
    void SetHiggsSignal(const bool higgsSignal);
    
    /// Set whether it is Drell-Yan sample
    void SetDrellYan(const bool isDrellYan);
    
    /// Set Drell-Yan decay channel for selection, access decay mode from nTuple branch
    void SetTrueLevelDYChannel(const int dy);
    
    /// Set name of output file
    void SetOutputfilename(const TString& outputfilename);
    
    /// Bool for separating direct dileptonic ttbar decays and decays via intermediate taus
    void SetRunViaTau(const bool runViaTau);
    
    /// Set folder for basic analysis output
    void SetAnalysisOutputBase(const char* analysisOutputBase);
    
    /// Set the Met recoil corrector
    void SetMetRecoilCorrector(const MetRecoilCorrector* const metRecoilCorrector);
    
    /// Whether MVA MET is used
    const bool& useMvaMet()const{return mvaMet_;}
    
    /// Set the kinematic reconstruction
    void SetKinematicReconstruction(KinematicReconstruction* kinematicReconstruction,
                                    const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors);
    
    /// Set the pileup reweighter
    void SetPileupScaleFactors(const PileupScaleFactors* pileupScaleFactors);
    
    /// Set the lepton scale factors
    void SetLeptonScaleFactors(const LeptonScaleFactors& scaleFactors);
    
    /// Set the trigger scale factors
    void SetTriggerScaleFactors(const TriggerScaleFactors& scaleFactors);
    
    /// Set the btag scale factors
    void SetBtagScaleFactors(BtagScaleFactors& scaleFactors);
    
    /// Set the MC sample which should be used for production of btag efficiencies
    void SetSampleForBtagEfficiencies(const bool isSampleForBtagEfficiencies);
    
    /// Set jet energy resolution scale factors
    void SetJetEnergyResolutionScaleFactors(const JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors);
    
    /// Set jet energy scale scale factors
    void SetJetEnergyScaleScaleFactors(const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors);
    
    /// Set top-pt scale factors
    void SetTopPtScaleFactors(const TopPtScaleFactors* topPtScaleFactors);
    
    /// Set histogram containing the number of weighted events in full sample
    void SetWeightedEvents(TH1* weightedEvents);

    /// Get version
    //no idea why this is needed! But Process() isn't called otherwise! Argh!
    virtual int Version() const { return 3; }

    /// Class definition
    ClassDef(AnalysisBase, 0);
    
    
    
protected:
    
    
    
// ----------------------- Protected methods for accessing the object structs -----------------------
    
    /// Get a constant reference to nTuple branches relevant for reconstruction level
    const RecoObjects& getRecoObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to nTuple branches holding generator information for all MC samples
    const CommonGenObjects& getCommonGenObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to nTuple branches for Top signal samples on generator level
    const TopGenObjects& getTopGenObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to nTuple branches for Higgs signal samples on generator level
    const HiggsGenObjects& getHiggsGenObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to nTuple branches for Z signal samples on generator level
    const ZGenObjects& getZGenObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to nTuple branches for kinematic reconstruction
    const KinRecoObjects& getKinRecoObjects(const Long64_t& entry)const;
    
    /// Get a constant reference to on-the-fly kinematic reconstruction results
    const KinRecoObjects& getKinRecoObjectsOnTheFly(const int leptonIndex, const int antiLeptonIndex,
                                                    const std::vector<int>& jetIndices, const std::vector<int>& bjetIndices,
                                                    const VLV& allLeptons, const VLV& jets, const std::vector<double>& jetBTagCSV,
                                                    const LV& met);
    
    
    /// Set for all object structs, that the specific nTuple entry is not read
    void resetObjectStructEntry()const;
    
    
    
// ----------------------- Protected methods for event and object selection -----------------------
    
    /// Check if opposite-charge dilepton combination exists,
    /// and check if lepton pair is correct flavour combination for the specified analysis channel (ee, emu, mumu)
    bool hasLeptonPair(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                       const std::vector<int>& lepPdgId)const;
    
    /// Check if event was triggered with the same dilepton trigger as the specified analysis channel
    bool failsDileptonTrigger(const Long64_t& entry)const;
    
    /// Set the b-tag algoritm and working point used for b-tag corrections and also for the discriminator cut value
    void setBtagAlgorithmAndWorkingPoint(const Btag::Algorithm& algorithm,
                                         const Btag::WorkingPoint& workingPoint);
    
    /// Returns the b-tag discriminator cut value associated to the algorithm and working point set in the b-tag tool
    double btagCutValue()const;
    
    /// Select events from Drell-Yan samples which need to be removed due to generator selection
    bool failsDrellYanGeneratorSelection(const Long64_t& entry)const;

    /// Select events from Top signal samples which need to be removed due to generator selection
    bool failsTopGeneratorSelection(const Long64_t& entry)const;
    
    /// Set usage of MVA MET instead of PF MET
    void mvaMet();
    
    
    
    
    
// ----------------------- Protected methods for application of corrections (e.g. scale factors) stored in the ntuple -----------------------
    
    /// Correct branching ratios of W decays in MadGraph samples
    double madgraphWDecayCorrection(const Long64_t& entry)const;
    
    /// Get weight due to pileup reweighting
    double weightPileup(const Long64_t& entry)const;
    
    /// Get weight due to generator weights
    double weightGenerator(const Long64_t& entry)const;
    
    /// Get weight of PDF variation
    double weightPdf(const Long64_t& entry, const int pdfNo)const;
    
    
    
// ----------------------- Protected methods for application of corrections (e.g. scale factors) NOT stored in the ntuple -----------------------
    
    ///  Recoil Correction of the MVA MET
    void correctMvaMet(const LV& dilepton, const int nJet, const Long64_t& entry)const;

    /// Get weight due to top-pt reweighting
    double weightTopPtReweighting(const Long64_t& entry)const;
    
    /// Get weight due to lepton efficiency MC-to-data scale factors
    double weightLeptonSF(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                          const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
    
    /// Get weight due to trigger efficiency MC-to-data scale factors
    double weightTriggerSF(const int leptonXIndex, const int leptonYIndex,
                           const VLV& allLeptons)const;
    
    /// Get weight due to b-tagging efficiency MC-to-data scale factors
    double weightBtagSF(const std::vector<int>& jetIndices,
                        const VLV& jets, const std::vector<int>& jetPartonFlavour,
                        const std::vector<double>& btagDiscriminators)const;
    
    /// Apply b-tag efficiency MC correction using random number based tag flipping
    void retagJets(std::vector<int>& bjetIndices, const std::vector<int>& jetIndices,
                   const VLV& jets, const std::vector<int>& jetPartonFlavours,
                   const std::vector<double>& btagDiscriminants)const;
    
    /// Get weight due to efficiency of kinematic reconstruction
    double weightKinReco()const;
    
    /// Whether to or not to produce b-tag efficiencies
    /// If it is true, no analysis output is produced, but only files containing the efficiencies
    bool makeBtagEfficiencies()const;
    
    /// Book histograms for b-tag efficiencies, in case the specified b-tag correction mode requires it
    void bookBtagEfficiencyHistos();
    
    /// Fill histograms for b-tag efficiencies, in case the specified b-tag correction mode requires it
    void fillBtagEfficiencyHistos(const std::vector<int>& jetIndices,
                                  const std::vector<double>& btagDiscriminators,
                                  const VLV& jets,
                                  const std::vector<int>& jetPartonFlavours,
                                  const double& weight);
    
    /// Produce b-tag efficiencies, in case the specified b-tag correction mode requires it
    void produceBtagEfficiencies();
    
    
    
    
// ----------------------- Various protected helper methods -----------------------
    
    /** Return a string describing the true level W+/W- decays from the ttbar system
     *
     * @return a string like e/tau->mu describing the decay to the W+/W- from the top/tbar decay
     *
     * Possible strings are:
     *   e/e  for the ee channel
     *   e/tau->mu for the W+ -> e+ decay and the W- -> tau- -> mu- decay
     *   ... and so on. The first part is from the top/W+, the second from the tbar/W-
     */
    const std::string topDecayModeString()const;

    /** Calculate the kinematic reconstruction and return whether at least one solution exists
     *
     * reconstruct the top quarks and store the result in the Hyp* member variables,
     * but also in the pointer kinRecoObjects
     *
     * @param leptonIndex Index of selected lepton in the vector of allLeptons
     * @param antiLeptonIndex Index of selected antilepton in the vector of allLeptons
     * @param jetIndices Indices of selected jets in the vector of jets
     * @param bjetIndices Indices of selected b-tagged jets in the vector of jets
     * @param allLeptons Vector of Lorentz vectors of all leptons from nTuple
     * @param jets Vector of Lorentz vectors of all jets from nTuple
     * @param jetBTagCSV B-tag discriminator values of all jets from nTuple
     * @param met Lorentz vector of MET
     *
     * @return true if there is at least one solution
     */
    bool calculateKinReco(const int leptonIndex, const int antiLeptonIndex, const std::vector<int>& jetIndices, const std::vector<int>& bjetIndices,
                          const VLV& allLeptons, const VLV& jets, const std::vector<double>& jetBTagCSV, const LV& met);
    
    /// Get H_t of jets
    double getJetHT(const std::vector<int>& jetIndices, const VLV& jets)const;
    
    /// Access identifier key of Higgs decay mode
    int higgsDecayMode(const Long64_t& entry)const;
    
    /// Store the object in the output list and return it
    template<class T> T* store(T* obj){return common::store(obj, fOutput);}
    
    /// Returns the samplename which is written as metadata in ingoing and outgoing root-file
    const TString& samplename()const{return samplename_;}
    
    /// Returns the output filename
    const TString& outputFilename()const{return outputfilename_;}
    
    /// Returns the product of the PDG IDs of the two leptons forming the dilepton system
    const int& channelPdgIdProduct()const{return channelPdgIdProduct_;}
    
    /// Returns the decay channel
    const Channel::Channel& channel()const{return channel_;}
    
    /// Returns the analysed systematic
    const Systematic::Systematic& systematic()const{return systematic_;}
    
    /// Whether it is a MC sample
    const bool& isMC()const{return isMC_;}
    
    /// Whether generator information about a ttbar system is stored (and whether it should be used for analysis)
    const bool& isTopSignal()const{return isTopSignal_;}
    
    /// Whether generator information about a Higgs system is stored (and whether it should be used for analysis)
    const bool& isHiggsSignal()const{return isHiggsSignal_;}
    
    /// Whether it is a Drell-Yan sample
    const bool& isDrellYan()const{return isDrellYan_;};
    
    /// Whether it is a ttbar sample (not ttbar + something)
    const bool& isTtbarSample()const{return isTtbarSample_;}
    
    /// Whether it is the ttbarsignalplustau sample (the part of ttbar which contains generator information about the ttbar system)
    const bool& isTtbarPlusTauSample()const{return isTtbarPlusTauSample_;}
    
    
    
    
    
private:
    
    
    
// ----------------------- Private methods for accessing the values in the nTuple -----------------------
    
    /// Access event entry for nTuple branches relevant for reconstruction level
    void GetRecoBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches holding generator information for all MC samples
    void GetCommonGenBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches of kinematic reconstruction
    void GetKinRecoBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches for Top signal samples on generator level
    void GetTopSignalBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches for Higgs signal samples on generator level
    void GetHiggsSignalBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches for Z signal samples on generator level
    void GetZSignalBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches for trigger bits
    void GetTriggerBranchesEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch of true vertex multiplicity
    void GetVertMultiTrueEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch for generator event weight
    void GetWeightGeneratorEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch for PDF weights
    void GetPDFEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch for Top decay mode
    void GetTopDecayModeEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch for Higgs decay mode
    void GetHiggsDecayModeEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branch for Z decay mode
    void GetZDecayModeEntry(const Long64_t& entry)const;
    
    /// Access event entry for nTuple branches of generated top and anti-top quark, used for reweighting
    void GetGenTopBranchesEntry(const Long64_t& entry)const;
    
    
    
// ----------------------- Private methods for setting the pointers to the nTuple -----------------------
    
    /// Clear all branches
    void clearBranches();
    
    /// Clear all variables assigned to branches
    void clearBranchVariables();
    
    
    /// Set addresses of nTuple branches relevant for reconstruction level
    void SetRecoBranchAddresses();
    
    /// Set address of nTuple branches for trigger bits
    void SetTriggerBranchAddresses();
    
    /// Set addresses of nTuple branches holding generator information for all MC samples
    void SetCommonGenBranchAddresses();
    
    /// Set addresses of nTuple branches of kinematic reconstruction
    void SetKinRecoBranchAddresses();
    
    /// Set address of nTuple branch of true vertex multiplicity
    void SetVertMultiTrueBranchAddress();
    
    /// Set address of nTuple branch for generator event weight
    void SetWeightGeneratorBranchAddress();
    
    /// Set address of nTuple branch for PDF weights
    void SetPdfBranchAddress();
    
    /// Set address of nTuple branch for Top decay mode
    void SetTopDecayBranchAddress();
    
    /// Set address of nTuple branch for Higgs decay mode
    void SetHiggsDecayBranchAddress();
    
    /// Set address of nTuple branch for Z decay mode
    void SetZDecayBranchAddress();
    
    /// Set addresses of nTuple branches for generated top and anti-top quark, used for reweighting
    void SetGenTopBranchAddresses();
    
    /// Set addresses of nTuple branches for Top signal samples on generator level
    void SetTopSignalBranchAddresses();
    
    /// Set addresses of nTuple branches for Higgs signal samples on generator level
    void SetHiggsSignalBranchAddresses();
    
    /// Set addresses of nTuple branches for Z signal samples on generator level
    void SetZSignalBranchAddresses();
    
    
    
    
// ----------------------- Private helper methods -----------------------
    
    /// Write analysis output to root file
    void writeOutput();
    
    
    
    
// ----------------------- Pointers to the nTuple and all its branches -----------------------
    
    /// Pointer to the analyzed TTree or TChain
    TTree* chain_;
    
    
    /// nTuple branches relevant for reconstruction level
    // Concerning physics objects
    TBranch* b_lepton;
    TBranch* b_lepPdgId;
    TBranch* b_lepID;
    TBranch* b_lepPfIso;
    TBranch* b_lepChargedHadronIso;
    TBranch* b_lepNeutralHadronIso;
    TBranch* b_lepPhotonIso;
    TBranch* b_lepPuChargedHadronIso;
    TBranch* b_lepCombIso;
    TBranch* b_lepDxyVertex0;
    TBranch* b_lepDzVertex0;
    TBranch* b_lepTrigger;
    TBranch* b_jet;
    TBranch* b_jetBTagTCHE;
    TBranch* b_jetBTagTCHP;
    TBranch* b_jetBTagSSVHE;
    TBranch* b_jetBTagSSVHP;
    TBranch* b_jetBTagJetProbability;
    TBranch* b_jetBTagJetBProbability;
    TBranch* b_jetBTagCSV;
    TBranch* b_jetBTagCSVMVA;
    TBranch* b_jetChargeGlobalPtWeighted;
    TBranch* b_jetChargeRelativePtWeighted;
    TBranch* b_jetPfCandidateTrack;
    TBranch* b_jetPfCandidateTrackCharge;
    TBranch* b_jetPfCandidateTrackId;
    TBranch* b_jetPfCandidateTrackIndex;
    TBranch* b_jetSelectedTrack;
    TBranch* b_jetSelectedTrackIPValue;
    TBranch* b_jetSelectedTrackIPSignificance;
    TBranch* b_jetSelectedTrackCharge;
    TBranch* b_jetSelectedTrackIndex;
    TBranch* b_jetSelectedTrackMatchToPfCandidateIndex;
    TBranch* b_jetSecondaryVertex;
    TBranch* b_jetSecondaryVertexPtCorrectedMass;
    TBranch* b_jetSecondaryVertexJetIndex;
    TBranch* b_jetSecondaryVertexFlightDistanceValue;
    TBranch* b_jetSecondaryVertexFlightDistanceSignificance;
    TBranch* b_jetSecondaryVertexTrackVertexIndex;
    TBranch* b_jetSecondaryVertexTrackMatchToSelectedTrackIndex;
    TBranch* b_met;
    TBranch* b_jetForMET;
    
    // Concerning event
    TBranch* b_runNumber;
    TBranch* b_lumiBlock;
    TBranch* b_eventNumber;
    
    
    /// nTuple branches holding trigger bits
    TBranch* b_triggerBits;
    TBranch* b_triggerBitsTau;
    TBranch* b_firedTriggers;
    
    
    /// nTuple branches holding generator information for all MC samples
    // Concerning physics objects
    TBranch* b_jetJERSF;
    TBranch* b_jetForMETJERSF;
    TBranch* b_vertMulti;
    TBranch* b_associatedGenJet;
    TBranch* b_associatedGenJetForMET;
    TBranch* b_jetPartonFlavour;
    TBranch* b_jetPartonFlavourForMET;
    
    
    /// nTuple branches of kinematic reconstruction
    TBranch* b_HypTop;
    TBranch* b_HypAntiTop;
    TBranch* b_HypLepton;
    TBranch* b_HypAntiLepton;
    TBranch* b_HypNeutrino;
    TBranch* b_HypAntiNeutrino;
    TBranch* b_HypB;
    TBranch* b_HypAntiB;
    TBranch* b_HypWPlus;
    TBranch* b_HypWMinus;
    TBranch* b_HypJet0index;
    TBranch* b_HypJet1index;
    
    
    /// nTuple branch of true vertex multiplicity
    TBranch* b_vertMultiTrue;
    
    
    /// nTuple branch for generator event weight
    TBranch* b_weightGenerator;
    
    
    /// nTuple branch for PDF weights
    TBranch* b_weightPDF;
    
    
    /// nTuple branch for Top decay mode
    TBranch* b_TopDecayMode;
    
    
    /// nTuple branch for Higgs decay mode
    TBranch* b_HiggsDecayMode;
    
    
    /// nTuple branch for Z decay mode
    TBranch* b_GenZDecayMode;
    
    
    /// nTuple branches for Top signal samples on generator level
    TBranch* b_GenTop;
    TBranch* b_GenAntiTop;
    TBranch* b_GenLepton;
    TBranch* b_GenAntiLepton;
    TBranch* b_GenLeptonPdgId;
    TBranch* b_GenAntiLeptonPdgId;
    TBranch* b_GenTau;
    TBranch* b_GenAntiTau;
    TBranch* b_GenNeutrino;
    TBranch* b_GenAntiNeutrino;
    TBranch* b_GenB;
    TBranch* b_GenAntiB;
    TBranch* b_GenWPlus;
    TBranch* b_GenWMinus;
    TBranch* b_GenMet;
    TBranch* b_allGenJets;
    TBranch* b_GenParticleP4;
    TBranch* b_GenParticlePdgId;
    TBranch* b_GenParticleStatus;
    TBranch* b_BHadJetIndex;
    TBranch* b_AntiBHadJetIndex;
    TBranch* b_BHadrons;
    TBranch* b_AntiBHadrons;
    TBranch* b_BHadronFromTopB;
    TBranch* b_AntiBHadronFromTopB;
    TBranch* b_BHadronVsJet;
    TBranch* b_AntiBHadronVsJet;
    TBranch* b_jetAssociatedPartonPdgId;
    TBranch* b_jetAssociatedParton;
    
    TBranch* b_genBHadPlusMothersPdgId;
    TBranch* b_genBHadPlusMothersStatus;
    TBranch* b_genBHadPlusMothersIndices;
    TBranch* b_genBHadPlusMothers;
    TBranch* b_genBHadIndex;
    TBranch* b_genBHadFlavour;
    TBranch* b_genBHadJetIndex;
    TBranch* b_genBHadLeptonIndex;
    TBranch* b_genBHadLeptonHadronIndex;
    TBranch* b_genBHadLeptonViaTau;
    TBranch* b_genBHadFromTopWeakDecay;
    
    TBranch* b_genCHadJetIndex;
    TBranch* b_genCHadLeptonIndex;
    TBranch* b_genCHadLeptonHadronIndex;
    TBranch* b_genCHadLeptonViaTau;
    TBranch* b_genCHadFromBHadron;
    
    TBranch* b_genExtraTopJetNumberId;
    
    
    /// nTuple branches for Higgs signal samples on generator level
    TBranch* b_GenH;
    TBranch* b_GenBFromH;
    TBranch* b_GenAntiBFromH;
    
    
    /// nTuple branches for Z signal samples on generator level
    TBranch* b_GenZ;
    TBranch* b_GenZMeDaughterParticle;
    TBranch* b_GenZMeDaughterAntiParticle;
    TBranch* b_GenZStableLepton;
    TBranch* b_GenZStableAntiLepton;
    
    
    
    
// ----------------------- Variables associated to the nTuple branches -----------------------
    
    /// Struct for holding variables associated to nTuple branches relevant for reconstruction level
    RecoObjects* recoObjects_;
    
    /// Struct for holding variables associated to nTuple branches holding generator information for all MC samples
    CommonGenObjects* commonGenObjects_;
    
    /// Struct for holding variables associated to nTuple branches for Top signal samples on generator level
    TopGenObjects* topGenObjects_;
    
    /// Struct for holding variables associated to nTuple branches for Higgs signal samples on generator level
    HiggsGenObjects* higgsGenObjects_;
    
    /// Struct for holding variables associated to nTuple branches for Z signal samples on generator level
    ZGenObjects* zGenObjects_;
    
    /// Struct for holding variables associated to nTuple branches of kinematic reconstruction
    KinRecoObjects* kinRecoObjects_;
    
    
    /// Variables associated to nTuple branches holding trigger bits
    UInt_t triggerBits_;
    //UInt_t triggerBitsTau_;
    //std::vector<std::string>* firedTriggers_;
    
    /// Variables associated to nTuple branch of true vertex multiplicity
    Int_t vertMultiTrue_;
    
    /// Variables associated to nTuple branch for generator event weight
    Double_t weightGenerator_;
    
    /// Variables associated to nTuple branch for PDF weights
    std::vector<double>* weightPDF_;
    
    /// Variables associated to nTuple branch for Top decay mode
    Int_t topDecayMode_;
    
    /// Variables associated to nTuple branch for Higgs decay mode
    Int_t higgsDecayMode_;
    
    /// Variables associated to nTuple branch for Z decay mode
    std::vector<int>* v_genZDecayMode_;
    
    
    
// ----------------------- Variables associated to nTuple metadata -----------------------
    
    /// Histogram holding the number of weighted events of the full sample as stored in nTuple, configured from outside
    TH1* h_weightedEvents;
    
    /// Information in nTuple stored in TObjString once per file, but added from outside and potentially configured
    TString samplename_;
    Channel::Channel channel_;
    Systematic::Systematic systematic_;
    bool isMC_;
    bool isTopSignal_;
    bool isHiggsSignal_;
    
    
    
// ----------------------- Helper variables needed for AnalysisBase -----------------------
    
    /// Further variables added from the outside
    bool isDrellYan_;
    bool isTtbarPlusTauSample_;
    bool correctMadgraphBR_;
    int channelPdgIdProduct_;
    #ifndef __CINT__
    std::function<bool(Long64_t)> checkZDecayMode_;
    #endif
    TString outputfilename_;
    bool runViaTau_;
    
    /// Whether it is a ttbar sample (and not ttbarH, ttbarW, ttbarZ, or any other thing)
    bool isTtbarSample_;
    
    /// Event counter
    Int_t eventCounter_;
    
    /// Directory where to store the analysis output
    const char* analysisOutputBase_;
    
    
    
    /// Pointer to the kinematic reconstruction instance
    KinematicReconstruction* kinematicReconstruction_;
    
    /// Pointer to the kinematic reconstruction scale factors instance
    const KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_;
    
    /// Pointer to the pileup reweighter instance
    const PileupScaleFactors* pileupScaleFactors_;
    
    /// Pointer to lepton scale factors instance
    const LeptonScaleFactors* leptonScaleFactors_;
    
    /// Pointer to trigger scale factors instance
    const TriggerScaleFactors* triggerScaleFactors_;
    
    /// Pointer to jet energy resolution scale factors instance
    const JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors_;
    
    /// Pointer to jet energy resolution scale factors instance
    const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors_;
    
    /// Pointer to btag scale factors instance
    BtagScaleFactors* btagScaleFactors_;
    
    /// Pointer to the top-pt scale factors instance
    const TopPtScaleFactors* topPtScaleFactors_;
    
    /// Whether the sample should be used for production of btag efficiencies
    bool isSampleForBtagEfficiencies_;
    
    /// Whether to use MVA MET
    bool mvaMet_;
    
    /// Pointer to the Met Recoil correction tools
    const MetRecoilCorrector* metRecoilCorrector_;
};





#endif




