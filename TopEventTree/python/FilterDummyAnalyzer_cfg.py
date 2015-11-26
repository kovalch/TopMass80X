import FWCore.ParameterSet.Config as cms

process = cms.Process("FilterDummyAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.FilterDummyAnalyzer = cms.EDAnalyzer('FilterDummyAnalyzer'
)


process.p = cms.Path(process.FilterDummyAnalyzer)
