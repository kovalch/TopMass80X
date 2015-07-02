import FWCore.ParameterSet.Config as cms

EventWeightMCSystematic = cms.EDProducer("EventWeightMCSystematic",
                                         genEventInfoTag=cms.InputTag("generator"),
                                         lheEventInfoTag=cms.InputTag("externalLHEProducer"),
                                         weightID=cms.string("0000"),
                                         printLHE=cms.bool(False)
                                         )


