import FWCore.ParameterSet.Config as cms

process = cms.Process("JetPropertiesProducerTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('JetPropertiesProducer')
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


#process.load("TopAnalysis.Configuration.Summer12.TTH_HToBB_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.load("TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.maxEvents.input = 11


process.load("TopAnalysis.HiggsUtils.producers.JetPropertiesProducer_cfi")
#process.jetProperties.src = ""



# Needs to be defined before PAT sequence
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_jetPropertiesProducer.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_*_*_JetPropertiesProducerTest',
    ),
)


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START53_V7G::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# Testing requires PAT objects
pfpostfix = ""
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *
jetCorr =('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=pfpostfix, jetCorrections=jetCorr, pvCollection=cms.InputTag('offlinePrimaryVertices'),typeIMetCorrections=True) 

from PhysicsTools.PatAlgos.tools.coreTools import *
# remove cleaning as it is not used (seems to be already applied in PF2PAT, so no need for this)
#removeCleaning(process, outputModules = [], postfix = pfpostfix)

# switch off MC matching
#removeMCMatching(process, ['All'])

# Remove some unneeded objects
#removeSpecificPATObjects( process, ['Taus', 'Photons', 'Muons', 'Electrons'] )
removeSpecificPATObjects( process, ['Taus', 'Photons', 'Muons', 'Electrons'], outputModules = [], postfix = pfpostfix )
#process.patDefaultSequence.remove( process.patTaus )


## remove the full pftau sequence as it is not needed for us
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauPFJets08Region'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauPileUpVertices'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauTagInfoProducer'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsPiZeros'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsLegacyTaNCPiZeros'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsLegacyHPSPiZeros'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBase'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsSelectionDiscriminator'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauProducerSansRefs'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauProducer'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBaseDiscriminationByDecayModeFinding'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBaseDiscriminationByLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTaus'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfNoTau'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByDecayModeFinding'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseChargedIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseChargedIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumChargedIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightChargedIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightIsolation'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseElectronRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumElectronRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightElectronRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMVAElectronRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseMuonRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumMuonRejection'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightMuonRejection'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFCandidates'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFChargedHadrons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFNeutralHadrons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFGammas'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patTaus'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatTaus'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatTaus'+pfpostfix))

## removal of unnecessary modules
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patCandidateSummary'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatCandidateSummary'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatElectrons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatMuons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatLeptons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatJets'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatPFParticles'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfPhotonSequence'+pfpostfix))

massSearchReplaceAnyInputTag(getattr(process,'patPF2PATSequence'+pfpostfix),'pfNoTau'+pfpostfix,'pfJets'+pfpostfix)






process.p = cms.Path(
    getattr(process,'patPF2PATSequence'+pfpostfix) *
    process.jetProperties)





process.e = cms.EndPath(process.out)
