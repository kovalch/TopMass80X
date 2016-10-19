#
# This file contains the Top PAG reference selection for electron + jets analysis.
#


### ------------------------- Reference selection -------------------------- ###


### Trigger selection

# HLT selection
triggerSelectionDataElectron  = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v* ' 
triggerSelectionMCElectron    = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v*' 

triggerSelectionDataMuon  = 'HLT_IsoMu18_v*' 
triggerSelectionMCMuon    = 'HLT_IsoMu18_v*' 

### Muon selection

# Minimal selection for all muons, also basis for signal and veto muons
# Muon ID ("loose")
muonCut  =     'isPFMuon'                                                                      # general reconstruction property
muonCut += ' && (isGlobalMuon || isTrackerMuon)'                                               # general reconstruction property
# Kinematics
muonCut += ' && pt > 15.'                                                                      # transverse momentum
muonCut += ' && abs(eta) < 2.4'                                                                # pseudo-rapisity range
# (Relative) isolation
#muonCut += ' && (chargedHadronIso+neutralHadronIso+photonIso-0.5*puChargedHadronIso)/pt < 0.2' # relative isolation w/ Delta beta corrections (factor 0.5)
muonCut += ' && (pfIsolationR04().sumChargedHadronPt + max( pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt,0.0)) / pt() <0.25'

# Signal muon selection on top of 'muonCut'
# Muon ID ("tight")
signalMuonCut  =     'isPFMuon'                                                                               # general reconstruction property
signalMuonCut += ' && isGlobalMuon'                                                                           # general reconstruction property
signalMuonCut += ' && globalTrack.normalizedChi2 < 10.'                                                       # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && track.hitPattern.trackerLayersWithMeasurement > 5'                                      # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && globalTrack.hitPattern.numberOfValidMuonHits > 0'                                       # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && abs(dB) < 0.2'                                                                          # 2-dim impact parameter with respect to beam spot (s. "PAT muon configuration" above)
signalMuonCut += ' && innerTrack.hitPattern.numberOfValidPixelHits > 0'                                       # tracker reconstruction
signalMuonCut += ' && numberOfMatchedStations > 1'                                                            # muon chamber reconstruction
# Kinematics
signalMuonCut += ' && pt > 25.'                                                                               # transverse momentum
signalMuonCut += ' && abs(eta) < 2.1'                                                                         # pseudo-rapisity range
# (Relative) isolation
#signalMuonCut += ' && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso))/pt < 0.12' # relative isolation w/ Delta beta corrections (factor 0.5)
signalMuonCut += ' && (pfIsolationR04().sumChargedHadronPt + max( pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt,0.0)) / pt() <0.15'

muonVertexMaxDZ = 0.5 # DeltaZ between muon vertex and PV

### Jet selection

# Signal jet selection
# Jet ID
jetCut  =     'numberOfDaughters > 1'                                 # PF jet ID:
jetCut += ' && neutralHadronEnergyFraction < 0.99'                    # PF jet ID:
jetCut += ' && neutralEmEnergyFraction < 0.99'                        # PF jet ID:
jetCut += ' && (chargedEmEnergyFraction < 0.99 || abs(eta) >= 2.4)'   # PF jet ID:  warum '|| abs(eta) >= 2.4)' ist das schneller?
jetCut += ' && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4)' # PF jet ID:
jetCut += ' && (chargedMultiplicity > 0 || abs(eta) >= 2.4)'          # PF jet ID:
# Kinematics
jetCut += ' && abs(eta) < 2.4'                                        # pseudo-rapisity range
# varying jet pt thresholds
veryLooseJetCut = 'pt > 30.' # transverse momentum (4 jets)
looseJetCut     = 'pt > 30.' # transverse momentum (3 jets)
tightJetCut     = 'pt > 30.' # transverse momentum (2 jets)
veryTightJetCut = 'pt > 30.' # transverse momentum (leading jet)


##TODO lot of ElectronCuts for data duplicated in analyzeTopHypotheses_cfg ... still?
### Electron selection
#Signalelektron
# Electron ID
electronGsfCut  =     '(electronID("{0}-standalone-loose")==1. || electronID("{0}-standalone-loose")==3. || electronID("{0}-standalone-loose")==5. || electronID("{0}-standalone-loose")==7.)'   #FIXME check electronID for 76x
# Kinematics
electronGsfCut += ' && ecalDrivenMomentum.pt > 20.'              #FIXME ecal driven needed?                                             # transverse energy
electronGsfCut += ' && abs(ecalDrivenMomentum.eta) < 2.5'                                                                               # pseudo-rapisity range
# (Relative) isolation  #FIXME updater with 'effective Areas?'  electron selectionString?
electronGsfCut += ' && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.15' # relative isolation with Delta beta corrections
# ... using re-calibrated (with regression energy) kinematics
electronCalibCut = electronGsfCut.replace( 'ecalDrivenMomentum.', '' ) 
electronCut = electronGsfCut

signalElectronCut = 'pt > 30 && abs(eta) < 2.1' 


#ElektronVeto
# Minimal selection for veto electrons
# ... using GsfElectron kinematics
# Electron ID
electronGsfVetoCut  =     'electronID("{0}-standalone-veto")'                                                  # electrons ID
# Kinematics
electronGsfVetoCut += ' && ecalDrivenMomentum.pt > 15.'                                                                                     # transverse energy
electronGsfVetoCut += ' && abs(ecalDrivenMomentum.eta) < 2.4'                                                                               # pseudo-rapisity range
# (Relative) isolation #FIXME updater with 'effective Areas?'  electron selectionString?
electronGsfVetoCut += ' && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.15' # relative isolation with Delta beta corrections

# ... using re-calibrated (with regression energy) kinematics
electronCalibVetoCut = electronGsfVetoCut.replace( 'ecalDrivenMomentum.', '' ) 
electronVetoCut = electronGsfVetoCut
### ------------------------------------------------------------------------ ###

#dileptonElecVeto
dileptonElectronVetoCut = 'electronID("{0}-standalone-loose")' #formely on veto, recommendation states loose
#dileptonElectronVetoCut +=' && ecalDrivenMomentum.pt > 20. && abs(ecalDrivenMomentum.eta) < 2.5'    
#dileptonElectronVetoCut +=' && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.2'


### Electron Conversion Rejection
conversionRejectionCut = '(electronID("{0}-standalone-loose")==4. || electronID("{0}-standalone-loose")==5. || electronID("{0}-standalone-loose")==6. || electronID("{0}-standalone-loose")==7.)'  #check Electron.h

#reminder ElectronID value map
#0: fails
#1: passes electron ID only
#2: passes electron Isolation only
#3: passes electron ID and Isolation only
#4: passes conversion rejection
#5: passes conversion rejection and ID
#6: passes conversion rejection and Isolation
#7: passes the whole selection




# Signal b-tagged jet selection
bTagCut = 'bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8' #TODO sure?


### Trigger matching
# Trigger object selection
triggerObjectSelectionData = 'type("TriggerMuon") && ( path("%s") )'%( triggerSelectionDataMuon )
triggerObjectSelectionMC   = 'type("TriggerMuon") && ( path("%s") )'%( triggerSelectionMCMuon )
