import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of jets
#
analyzeJets = cms.EDAnalyzer("JetEventAnalyzer",
    jets            = cms.InputTag("goodJetsPF30"),
    alternativeJets = cms.InputTag(""),
    met = cms.InputTag("patMETs"),
    gluonTagSrc  = cms.InputTag(""),
    maxNJets = cms.int32(20)
)
