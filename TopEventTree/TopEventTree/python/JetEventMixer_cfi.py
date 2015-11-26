import FWCore.ParameterSet.Config as cms

#
# 
#

mixJets = cms.EDProducer("JetEventMixer",
                         input = cms.SecSource("PoolSource",
                                               fileNames =  cms.untracked.vstring()
                                               ),
                         ## number of events to be mixed
                         nMix    = cms.int32(10),
                         ## min number of jets to be present, which will be mixed
                         nMixMin = cms.int32(6),
                         ## speed up factor used to randomly selected 1:speedUp_ events
                         speedUp = cms.int32(1)
)
