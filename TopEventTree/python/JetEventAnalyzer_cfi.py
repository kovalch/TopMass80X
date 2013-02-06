import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of jets
#
analyzeJets = cms.EDAnalyzer("JetEventAnalyzer",
    jets         = cms.InputTag("goodJetsPF30"),
    #allJets      = cms.InputTag("patJets"),
    #noPtEtaJets  = cms.InputTag("noPtEtaJetsPF"),
    maxNJets = cms.int32(10)
)
