import FWCore.ParameterSet.Config as cms

filterWNJets = cms.EDFilter("GeneratorWNJetsFilter",
    NJet       = cms.int32(0),
)
