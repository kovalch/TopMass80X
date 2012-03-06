import FWCore.ParameterSet.Config as cms
import os
import sys

process = cms.Process("topMass")

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
       '/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0004/BA357784-557D-E011-A636-0017A477000C.root'
#       '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/832/22EBE93E-B1DB-E011-9C16-BCAEC518FF8A.root',
#       '/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/00454769-577B-E011-ACCD-001E0B49808A.root',
#       '/store/data/Run2011A/SingleMu/AOD/PromptReco-v4/000/166/462/8A8DCE75-8190-E011-B082-001D09F2915A.root'
#       '/store/data/Run2011A/SingleMu/AOD/PromptReco-v6/000/172/620/24054E7E-17C0-E011-AA64-001D09F28D4A.root'
       ] );



secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

## configure geometry & conditions
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')

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

## selection
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu17_v*","HLT_IsoMu24_v*","HLT_IsoMu24_eta2p1_v*"], throw=False)

process.leadingJetSelection.src = 'tightLeadingPFJets'
process.bottomJetSelection.src  = 'tightBottomPFJets'

## b-tag selection
process.tightBottomPFJets.cut = 'bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74';
process.bottomJetSelection.minNumber = 0;

from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "patMETsPF")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "simpleSecondaryVertexHighEffBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 1.74
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 1.74
process.hitFitTtSemiLepEventHypothesis.useBTagging       = False

## choose which hypotheses to produce
addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
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
                        process.hltFilter *
                        process.semiLeptonicSelection *
                        process.semiLeptonicEvents *
                        process.makeTtSemiLepEvent *
                        process.analyzeHypotheses
                        )

process.path.remove(process.centralJets)
process.path.remove(process.reliableJets)
process.path.remove(process.goodJets)
process.path.remove(process.trackCountingHighPurBJets)
process.path.remove(process.trackCountingHighEffBJets)
process.path.remove(process.tightLeadingJets)
process.path.remove(process.tightBottomJets)
process.path.remove(process.unconvTightElectronsEJ)
process.path.remove(process.goodElectronsEJ)
process.path.remove(process.looseElectronsEJ)
process.path.remove(process.tightElectronsEJ)
                        
from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
prependPF2PATSequence(process, options = {'runOnOLDcfg': True,
                                          'runOnMC': False,
                                          'runOnAOD': True,
                                          'electronIDs': '',
                                          'switchOffEmbedding': False,
                                          'pfIsoConeMuon': 0.4,
                                          'pfIsoConeElec': 0.4,
                                          'skipIfNoPFMuon': True,
                                          'METCorrectionLevel': 2,
                                          })

## adaptions (re-aranging of modules) to speed up processing
pathnames = process.paths_().keys()
for pathname in pathnames:
    ## move the trigger to the beginning of the sequence
    getattr(process, pathname).remove(process.hltFilter)
    getattr(process, pathname).insert(0,process.hltFilter)


