import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.categories.append('TtSemiLeptonicEvent')
process.MessageLogger.cerr.TtSemiLeptonicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

## define input
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/snaumann/TTJets_TuneD6T_7TeV-madgraph-tauola/SingleMuSkim_1/77c3ad1ea85046d30b430315aa7f3138/PAT_9_1_CJg.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_9_1_dqi.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_8_3_1Py.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_7_1_6CQ.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_6_3_JEo.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_5_1_6so.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_4_1_mQd.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_3_1_RMF.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_2_3_7oC.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_1_3_3cA.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_13_1_0Y9.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_12_1_tlA.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_11_1_ySV.root',
       '/store/user/snaumann/Mu/SingleMuSkim-Run2010B-Nov4ReReco_Mu15_2/1bfce60840ad652d060bf69d2d1ff328/PAT_10_3_CGj.root'
       ] );



secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry & conditions
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_38Y_V14::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")

## redefine veto jets to be sure it is also replaced when running on PF
from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJets
process.vetoJets.src="goodJetsPF30"
process.vetoJets.cut=''

## sequences for TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "patMETsPF")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)

process.TtSemiLepJetCombMVAFileSource = cms.ESSource("TtSemiLepJetCombMVAFileSource",
  ttSemiLepJetCombMVA = cms.FileInPath('TopMass/Configuration/data/TtSemiLepJetComb.mva')
)

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
#process.ttSemiLepJetPartonMatch.maxDist   = 0.3
process.ttSemiLepJetPartonMatch.maxNJets = -1

process.kinFitTtSemiLepEventHypothesis.useBTagging = False
# 1: Whad-mass, 2: Wlep-mass, 3: thad-mass, 4: tlep-mass, 5: nu-mass, 6: equal t-masses
process.kinFitTtSemiLepEventHypothesis.constraints = 1, 2, 6

findTtSemiLepJetCombMVA.maxNComb = 1;

## choose which hypotheses to produce
addTtSemiLepHypotheses(process,
                       ["kKinFit", "kHitFit", "kMVADisc"]
                       )
removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
process.load("TopMass.Analyzer.EventHypothesisAnalyzer_data_cff")

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTop.root')
)

## end path   
process.path = cms.Path(#process.patDefaultSequence *
                        process.semiLeptonicSelection *
                        process.semiLeptonicEvents *
                        process.makeTtSemiLepEvent *
                        process.analyzeHypotheses
                        )
