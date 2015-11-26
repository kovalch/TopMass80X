import FWCore.ParameterSet.Config as cms

FilterDummyAnalyzer = cms.EDAnalyzer('FilterDummyAnalyzer',
	CutFlagsBranchName= cms.string("CutFlags.")
)
