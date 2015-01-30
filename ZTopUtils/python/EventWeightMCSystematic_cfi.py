import FWCore.ParameterSet.Config as cms

EventWeightMCSystematic = cms.EDProducer("EventWeightMCSystematic",
                                         weightId=cms.string("0")
                                         )


