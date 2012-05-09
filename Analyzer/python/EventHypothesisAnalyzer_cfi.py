import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of top event hypothese
#
analyzeHypothesis = cms.EDAnalyzer("EventHypothesisAnalyzer",
    semiLepEvent = cms.InputTag("ttSemiLepEvent"),
    hypoClassKey = cms.InputTag("ttSemiLepHypMaxSumPtWMass","Key"),
    
    jets         = cms.InputTag("goodJetsPF30"),
    noPtEtaJets  = cms.InputTag("noPtEtaJetsPF"),
    leps         = cms.InputTag("tightMuons"),
    
    PUSrc        = cms.InputTag("addPileupInfo"),
    
    VertexSrc    = cms.InputTag("goodOfflinePrimaryVertices"),
    
    PUWeightSrc  = cms.InputTag("eventWeightPUsysNo", "eventWeightPU"),
    PUWeightUpSrc  = cms.InputTag("eventWeightPUsysUp", "eventWeightPUUp"),
    PUWeightDownSrc  = cms.InputTag("eventWeightPUsysDown", "eventWeightPUDown"),
    
    bWeightSrc   = cms.InputTag("bTagSFEventWeight"),
    bWeightSrc_bTagSFUp     = cms.InputTag("bTagSFEventWeightBTagSFUp"),
    bWeightSrc_bTagSFDown   = cms.InputTag("bTagSFEventWeightBTagSFDown"),
    bWeightSrc_misTagSFUp   = cms.InputTag("bTagSFEventWeightMisTagSFUp"),
    bWeightSrc_misTagSFDown = cms.InputTag("bTagSFEventWeightMisTagSFDown"),
    
    muWeightSrc  = cms.InputTag("effSFMuonEventWeight"),
    
    savePDFWeights = cms.bool(False),
)


