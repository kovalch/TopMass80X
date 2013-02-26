import FWCore.ParameterSet.Config as cms

#
# 
#

mixJets = cms.EDProducer("JetEventMixer",
                         input = cms.SecSource("PoolSource",
                                               fileNames =  cms.untracked.vstring()
                                               ),
                         nMix = cms.int32(10)
)
