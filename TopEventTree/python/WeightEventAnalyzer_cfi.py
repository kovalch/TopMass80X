import FWCore.ParameterSet.Config as cms

#
# module to make event weights for the analysis
#
analyzeWeights = cms.EDAnalyzer("WeightEventAnalyzer",
    mcWeight  = cms.double(1.0),

    puSrc     = cms.InputTag("addPileupInfo"),
    vertexSrc = cms.InputTag("goodOfflinePrimaryVertices"),
    
    puWeightSrc  = cms.InputTag("eventWeightPUsysNo", "eventWeightPU"),
    puWeightUpSrc  = cms.InputTag("eventWeightPUsysUp", "eventWeightPUUp"),
    puWeightDownSrc  = cms.InputTag("eventWeightPUsysDown", "eventWeightPUDown"),
    
    jets                    = cms.InputTag("goodJetsPF30"),
    bWeightSrc              = cms.InputTag("bTagSFEventWeight"),
    bWeightSrc_bTagSFUp     = cms.InputTag("bTagSFEventWeightBTagSFUp"),
    bWeightSrc_bTagSFDown   = cms.InputTag("bTagSFEventWeightBTagSFDown"),
    bWeightSrc_misTagSFUp   = cms.InputTag("bTagSFEventWeightMisTagSFUp"),
    bWeightSrc_misTagSFDown = cms.InputTag("bTagSFEventWeightMisTagSFDown"),
    
    triggerWeightSrc  = cms.InputTag("effSFMuonEventWeight"),
    
    muWeightSrc  = cms.InputTag("effSFMuonEventWeight"),
    elWeightSrc  = cms.InputTag("effSFElectronEventWeight"),
    
    genEventSrc    = cms.InputTag("generator"),
    ttEvent        = cms.InputTag("ttSemiLepEvent"),
    savePDFWeights = cms.bool(False),
    brCorrection   = cms.bool(False),
)
