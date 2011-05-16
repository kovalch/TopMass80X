import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of top event hypothese
#
analyzeHypothesis = cms.EDAnalyzer("EventHypothesisAnalyzer",
    semiLepEvent = cms.InputTag("ttSemiLepEvent"),
    hypoClassKey = cms.InputTag("ttSemiLepHypMaxSumPtWMass","Key"),
    jets         = cms.InputTag("goodJetsPF30"),
)


