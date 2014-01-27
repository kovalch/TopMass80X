import FWCore.ParameterSet.Config as cms


correctMuonEnergy = cms.EDProducer("muonRochesterCorrector",
                                   ### the type will be input and output type
                                   ## possible values are:
                                       # patMuons
                                       # recoMuons
                                       # pfMuons
                                   muonType = cms.string("pfMuons"),
                                   muonSrc = cms.InputTag("pfSelectedMuons"),
                                   isMC = cms.bool(True),
                                   randomSeed= cms.int32(12345),
                                   debug=cms.bool(False)
                                   )
