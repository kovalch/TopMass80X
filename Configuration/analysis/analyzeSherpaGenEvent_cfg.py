import FWCore.ParameterSet.Config as cms
import glob

process = cms.Process("analyzeSherpaGenEvent")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

inputfiles = []
hepmcfiles = glob.glob('*.hepmc2g')
for hepmcfile in hepmcfiles:
  inputfiles.append('file:' + hepmcfile)
print inputfiles

process.source = cms.Source("MCFileSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(inputfiles)
)

## define input
#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,
#                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
#)
#
#readFiles.extend( [
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1312_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1290_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_404_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_402_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1288_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1296_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_590_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1292_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_390_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1300_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_1298_fastreco_FASTSIM_HLT_PU.root',
#          '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_610_fastreco_FASTSIM_HLT_PU.root',
#          
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1325_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1389_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1423_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1465_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1349_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1435_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1463_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1359_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1379_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1409_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1407_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1473_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1441_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1391_fastreco_FASTSIM_HLT_PU.root',
#          #'/store/user/mseidel/TT_lund_8TeV-sherpa/job_1387_fastreco_FASTSIM_HLT_PU.root',
#          
#          #'/store/mc/Summer12_DR53X/WJetsToLNu_scaleup_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0046B53F-45F1-E111-9774-00266CF25C88.root'
#  ] )
#
#
#secFiles.extend( [
#               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.genParticles.src = 'source'

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("TopMass.TopEventTree.SherpaGenEventAnalyzer_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
   maxEventsToPrint = cms.untracked.int32(10),
   printVertex = cms.untracked.bool(False),
   src = cms.InputTag("genParticles")
)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeSherpaGenEvent.root')
)

# register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = ""

## end path   
process.path = cms.Path(
                        process.genParticles *
                        process.genJetParticles *
                        process.printTree *
                        process.ak5GenJets *
                        process.analyzeSherpaGenEvent
                        )


