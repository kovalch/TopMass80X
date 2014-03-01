import FWCore.ParameterSet.Config as cms

process = cms.Process("ReadJetCorrectorDB")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'
#process.GlobalTag.globaltag = 'START42_V17::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.source = cms.Source('EmptySource')

process.reader = cms.EDAnalyzer('JetCorrectorDBReader',
                                payloadName = cms.untracked.string('AK5PF'),
                                globalTag = cms.untracked.string('FT_53_V21_AN5'),
                                printScreen = cms.untracked.bool(False),
                                createTextFile = cms.untracked.bool(True)
                                )

process.p = cms.Path(process.reader)

