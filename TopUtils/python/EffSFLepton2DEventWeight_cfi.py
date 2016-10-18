import FWCore.ParameterSet.Config as cms

effSFLepton2DEventWeight = cms.EDProducer("EffSFLepton2DEventWeight",
  particles  = cms.InputTag(""), ## jet collection (after jet selection, before b-tagging)
  sysVar   = cms.string(""),                  ## "noSys", "combinedEffSFNormUpStat/Down", "combinedEffSFShapeUpEta(Pt)/Down", 
                                              ## "combinedEffSFNormUpSys/Down"
                                              ## "PUup", "PUdown"
                                              ## "flatTriggerSF"
  verbose  = cms.int32(  0),                  ## set to 1 if terminal text output is desired
  filename = cms.FileInPath("TopAnalysis/Configuration/data/MuonEffSF2D2012.root"),
  histoname = cms.string(" "),
  additionalSystErr = cms.double(0.),
  etaMaxVal = cms.double(5.),
)
