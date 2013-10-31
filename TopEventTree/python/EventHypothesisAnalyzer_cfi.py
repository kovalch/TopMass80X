import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of top event hypothese
#
analyzeHypothesis = cms.EDAnalyzer("EventHypothesisAnalyzer",
    ttEvent      = cms.InputTag("ttSemiLepEvent"),
    hypoClassKey = cms.InputTag("ttSemiLepHypMaxSumPtWMass","Key"),
    ttEventGen2  = cms.InputTag(""),
    jets         = cms.InputTag("goodJetsPF30"), # needed in fullHad channel for reco masses
    topBranchName= cms.string("top."),
    
    lepton       = cms.int32(13), # (positive) pdgId
   
    #leps         = cms.InputTag("tightMuons"),
    #mets         = cms.InputTag("scaledMET:scaledMETs"),
    #
    #maxNJets = cms.int32(10),
    maxCombo = cms.int32(10000)
)
