import FWCore.ParameterSet.Config as cms

#
# module to make event weights for the analysis
#
analyzeWeights = cms.EDAnalyzer("WeightEventAnalyzer",
    mcWeightSrc  = cms.InputTag("eventWeightMC"),

    puSrc        = cms.InputTag("addPileupInfo"),
    vertexSrc    = cms.InputTag("goodOfflinePrimaryVertices"),
    
    puWeightSrc  = cms.InputTag("eventWeightPUsysNo", "eventWeightPU"),
    puWeightUpSrc  = cms.InputTag("eventWeightPUsysUp", "eventWeightPUUp"),
    puWeightDownSrc  = cms.InputTag("eventWeightPUsysDown", "eventWeightPUDown"),
    
    jets                    = cms.InputTag("goodJetsPF30"),
    bWeightSrc              = cms.InputTag("bTagSFEventWeight"),
    bWeightSrc_bTagSFUp     = cms.InputTag("bTagSFEventWeightBTagSFUp"),
    bWeightSrc_bTagSFDown   = cms.InputTag("bTagSFEventWeightBTagSFDown"),
    bWeightSrc_misTagSFUp   = cms.InputTag("bTagSFEventWeightMisTagSFUp"),
    bWeightSrc_misTagSFDown = cms.InputTag("bTagSFEventWeightMisTagSFDown"),
    
    triggerWeightSrc  = cms.InputTag(""),
    
    muWeightSrc  = cms.InputTag("effSFMuonEventWeight"),
    elWeightSrc  = cms.InputTag("effSFElectronEventWeight"),
    
    genEventSrc    = cms.InputTag("generator"),
    savePDFWeights = cms.bool(False),
    
)
