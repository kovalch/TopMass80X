import FWCore.ParameterSet.Config as cms

JetDebug = cms.EDProducer('JetDebug',
    jets  = cms.InputTag("slimmedJets"),
    evtSolLabel  = cms.InputTag("")
)
