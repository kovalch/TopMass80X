import FWCore.ParameterSet.Config as cms

process = cms.Process("USERTEST")

# ================
#  Message Logger
# ================

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ==========================
#  Configuration of Process
# ==========================

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service("TFileService", fileName = cms.string("MC_PUDist.root") )

#### Input Files

#process.load("TopAnalysis.Configuration.samples.Spring11_WJets_PATtuple_cff")

#process.load("TopAnalysis.Configuration.samples.Summer11_TTJets_TuneZ2_7TeV_madgraph_tauola_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_WJetsToLNu_TuneZ2_7TeV_madgraph_tauola_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_20to30_BCtoE_TuneZ2_7TeV_pythia6_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_20to30_EMEnriched_TuneZ2_7TeV_pythia6_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_30to80_BCtoE_TuneZ2_7TeV_pythia6_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_30to80_EMEnriched_TuneZ2_7TeV_pythia_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_80to170_BCtoE_TuneZ2_7TeV_pythia_cff")
#process.load("TopAnalysis.Configuration.samples.Summer11_QCD_Pt_80to170_EMEnriched_TuneZ2_7TeV_pythia6_cff")
#process.load("TopAnalysis.Configuration.zprime_M500GeV_W5000MeV_Madgraph_Summer11_AOD_cff")
#process.load("TopAnalysis.Configuration.zprime_M750GeV_W7500MeV_Madgraph_Summer11_AOD_cff")
process.load("TopAnalysis.Configuration.Spring16_TT_TuneCUETP8M1_13TeV-powheg-pythia8_v2_cff")
#process.load("TopAnalysis.Configuration.samples.Fall11_TTJets_TuneZ2_7TeV_powheg_tauola_cff")
#process.load("TopAnalysis.Configuration.samples.Fall11_TTJets_TuneZ2_7TeV_mcatnlo_cff")

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())

#### Number of Events

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

#### Load default configuration

process.load("TopAnalysis.TopAnalyzer.MCPileUp_cfi")

#### Define path

process.p = cms.Path(process.MCPUDistribution)
