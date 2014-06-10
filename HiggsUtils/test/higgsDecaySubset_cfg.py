import FWCore.ParameterSet.Config as cms


process = cms.Process("HiggsDecayTest")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


#process.load("TopAnalysis.Configuration.Summer12.TTH_HToBB_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.load("TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.maxEvents.input = 101


process.load("TopAnalysis.HiggsUtils.producers.HiggsDecaySubset_cfi")
process.decaySubsetHiggs.fillMode = 'kME'
process.decaySubsetHiggs.addRadiation = False


process.p = cms.Path(process.decaySubsetHiggs)





process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_higgsDecaySubset.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
	'keep recoGenParticles_*_*_SIM',
	'keep recoGenParticles_*_*_HiggsDecayTest',
    ),
)

process.e = cms.EndPath(process.out)
