import FWCore.ParameterSet.Config as cms

process = cms.Process("signalProcessID")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## define input
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,
                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
)

readFiles.extend( [
          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1312_fastreco_FASTSIM_HLT_PU.root'
          #'/store/mc/Summer12_DR53X/WJetsToLNu_scaleup_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0046B53F-45F1-E111-9774-00266CF25C88.root'
  ] )


secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

## configure process options
#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
   maxEventsToPrint = cms.untracked.int32(10),
   printVertex = cms.untracked.bool(False),
   src = cms.InputTag("genParticles")
)

## end path   
process.path = cms.Path(
                        process.printTree
                        )


