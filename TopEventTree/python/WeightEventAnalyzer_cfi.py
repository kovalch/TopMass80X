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
    bWeightSrc_bTagCjetSFUp     = cms.InputTag("bTagSFEventWeightBTagCjetSFUp"),
    bWeightSrc_bTagCjetSFDown   = cms.InputTag("bTagSFEventWeightBTagCjetSFDown"),
    bWeightSrc_misTagSFUp   = cms.InputTag("bTagSFEventWeightMisTagSFUp"),
    bWeightSrc_misTagSFDown = cms.InputTag("bTagSFEventWeightMisTagSFDown"),
    
    triggerWeightSrc        = cms.InputTag("effSFMuonEventWeight"),
    triggerWeightSrcUp      = cms.InputTag("effSFMuonEventWeightUp"),
    triggerWeightSrcDown    = cms.InputTag("effSFMuonEventWeightDown"),
    
    bJESSrc_fNuUp     = cms.InputTag("bJESEventWeightFNuUp"),
    bJESSrc_fNuDown   = cms.InputTag("bJESEventWeightFNuDown"),
    bJESSrc_frag      = cms.InputTag("bJESEventWeightFrag"),
    bJESSrc_fragHard  = cms.InputTag("bJESEventWeightFragHard"),
    bJESSrc_fragSoft  = cms.InputTag("bJESEventWeightFragSoft"),
    
    muWeightSrc  = cms.InputTag("effSFMuonEventWeight"),
    elWeightSrc  = cms.InputTag("effSFElectronEventWeight"),
    
    genEventSrc    = cms.InputTag("generator"),
    lheEventSrc    = cms.InputTag("externalLHEProducer"),
    ttEvent        = cms.InputTag("ttSemiLepEvent"),
    savePDFWeights = cms.bool(False),
    brCorrection   = cms.bool(False),
    showLHEweightTypes = cms.bool(False),
    
    topPtSFa =  cms.double(0.159),
    topPtSFb =  cms.double(-0.00141),
    topPtSFthreshold =   cms.double(400)
    
)
 
