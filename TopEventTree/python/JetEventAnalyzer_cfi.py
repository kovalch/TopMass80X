import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of jets
#
analyzeJets = cms.EDAnalyzer("JetEventAnalyzer",
    jets            = cms.InputTag("goodJetsPF30"),
    alternativeJets = cms.InputTag(""),
    #allJets      = cms.InputTag("patJets"),
    #noPtEtaJets  = cms.InputTag("noPtEtaJetsPF"),
    gluonTagSrc  = cms.InputTag(""),
    maxNJets = cms.int32(20)
)
