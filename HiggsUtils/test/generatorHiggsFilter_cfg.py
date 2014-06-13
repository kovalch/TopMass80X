import FWCore.ParameterSet.Config as cms

process = cms.Process("GeneratorHiggsFilterTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('GeneratorHiggsFilter')
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


#process.load("TopAnalysis.Configuration.Summer12.TTH_HToBB_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.load("TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.maxEvents.input = 101


process.load("TopAnalysis.HiggsUtils.sequences.higgsGenEvent_cff")
#process.decaySubsetHiggs.fillMode = 'kME'
#process.decaySubsetHiggs.addRadiation = False

process.load("TopAnalysis.HiggsUtils.filters.GeneratorHiggsFilter_cfi")
process.generatorHiggsFilter.channels = ["b"]
#process.generatorHiggsFilter.channels = ["W", "tau"]
#process.generatorHiggsFilter.channels = ["none"]
#process.generatorHiggsFilter.channels = []
#process.generatorHiggsFilter.invert_selection = True



process.p = cms.Path(process.makeGenEvtHiggs*process.generatorHiggsFilter)




process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_generatorHiggsFilter.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoGenParticles_*_*_SIM',
        #'keep recoGenParticles_*_*_HiggsGenTest',
	#'keep HiggsGenEvent_*_*_HiggsGenTest',
	'keep *_*_*_GeneratorHiggsFilterTest',
    ),

)

process.e = cms.EndPath(process.out)
