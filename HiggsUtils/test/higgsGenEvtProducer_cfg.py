import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsGenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('HiggsGenEvent')
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


#process.load("TopAnalysis.Configuration.Summer12.TTH_HToBB_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.load("TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff")
process.maxEvents.input = 11


process.load("TopAnalysis.HiggsUtils.producers.HiggsDecaySubset_cfi")
process.decaySubsetHiggs.fillMode = 'kME'
process.decaySubsetHiggs.addRadiation = False
process.load("TopAnalysis.HiggsUtils.producers.HiggsGenEvtProducer_cfi")


process.p = cms.Path(process.decaySubsetHiggs*process.genEvtHiggs)




process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_higgsGenEvtProducer.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoGenParticles_*_*_SIM',
        #'keep recoGenParticles_*_*_HiggsGenTest',
        #'keep HiggsGenEvent_*_*_HiggsGenTest',
        'keep *_*_*_HiggsGenTest',
    ),

)

process.e = cms.EndPath(process.out)
