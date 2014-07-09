import FWCore.ParameterSet.Config as cms

#
# module to make simple gen analyses with Sherpa
#
analyzeSherpaGenEvent = cms.EDAnalyzer("SherpaGenEventAnalyzer",
    genJets = cms.InputTag("ak5GenJets"),
)
