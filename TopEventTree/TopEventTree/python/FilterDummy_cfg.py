import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.myProducerLabel = cms.EDProducer('FilterDummy'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('FilterDummyOut.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
