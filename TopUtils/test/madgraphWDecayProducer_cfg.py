import FWCore.ParameterSet.Config as cms


process = cms.Process("MadgraphWDecayProducerTest")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('MadgraphWDecayProducer')
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


inputScript = "TopAnalysis.Configuration.Summer12.TTJets_MSDecays_central_TuneZ2star_8TeV_madgraph_tauola__Summer12_DR53X_PU_S10_START53_V19_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WZJetsTo3LNu_TuneZ2_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WWJetsTo2L2Nu_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WWZNoGstarJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.TTWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
#inputScript = "TopAnalysis.Configuration.Summer12.WWGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
###inputScript = "TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
process.load(inputScript)
#process.source.skipEvents = cms.untracked.uint32(200000)
process.maxEvents.input = 101


if inputScript.find("_madgraph_") == -1:
    print 'MadGraph sample? --> NO'
    process.madgraphWDecaySequence = cms.Sequence()
else:
    print 'MadGraph sample? --> YES'
    process.load("TopAnalysis.TopUtils.MadgraphWDecayProducer_cfi")
    process.madgraphWDecaySequence = cms.Sequence(process.MadgraphWDecayProducer)


process.p = cms.Path(process.madgraphWDecaySequence)





process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_madgraphWDecayProducer.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        #'keep recoGenParticles_*_*_SIM',
        'keep *_*_*_MadgraphWDecayProducerTest',
    ),
)

process.e = cms.EndPath(process.out)
