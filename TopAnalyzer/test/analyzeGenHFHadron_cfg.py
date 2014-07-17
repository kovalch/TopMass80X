import FWCore.ParameterSet.Config as cms

## ---
##   use this file to test the GenHFHadronAnalyzer.cc module
## ---


# set sequence shortcut
process = cms.Process("Analyzer")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## define input
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(    
    ## add your favourite file here
    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00000/FC1E7E5A-E60F-E211-9605-485B39800C3B.root',
    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00000/FAB3BB25-2210-E211-968B-20CF305616E0.root',
    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00000/F8FDB54E-1110-E211-A162-E0CB4E19F962.root'
    ),
    skipEvents = cms.untracked.uint32(0)
 )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeGenHFHadronJetsAnalyzer_test.root')
)

## ---
##    load GenParticle
## ---

# get gen ttbar event
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

# get genJets with hadrons in it
process.load("TopAnalysis.TopUtils.sequences.improvedJetHadronQuarkMatching_cff")

# supply PDG ID to real name resolution of MC particles, necessary for GenHFHadronAnalyzer
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

## get particle content of sample with IDs
process.load("TopAnalysis.TopAnalyzer.GenHFHadronAnalyzer_cff")
process.analyzeGenHFHadronJets = process.analyzeGenHFHadronJets.clone()

process.analyzeGenHFHadronJets.noBBbarResonances = True


## ---
##    run the final sequence
## ---

process.p1 = cms.Path(
    ## apply the analyzer
    process.makeGenEvt *
    process.improvedJetHadronQuarkMatchingSequence *
    process.analyzeGenHFHadronJets
    )
