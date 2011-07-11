import FWCore.ParameterSet.Config as cms

lJesFactor = '@lJesFactor@'
if lJesFactor.startswith('@'):
  lJesFactor = '1.0'

bJesFactor = '@bJesFactor@'
if bJesFactor.startswith('@'):
  bJesFactor = '1.0'

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
       '/store/user/eschliec/TTJets_TuneD6T_7TeV-madgraph-tauola/PATWithPF_v4/e59efddd8a1547799dca5b47d5556447/patTuple_39_1_4mM.root'
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
process.GlobalTag.globaltag = cms.string('START41_V0::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")

## redefine veto jets to be sure it is also replaced when running on PF
from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJets
process.vetoJets.src="goodJetsPF30"
process.vetoJets.cut=''

## sequences for ttGenEvent and TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "leps", "tightMuons")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "patMETsPF")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"

## choose which hypotheses to produce
addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )
#removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
process.load("TopMass.Analyzer.EventHypothesisAnalyzer_cff")

## configure JetEnergyScale
process.hitFitTtSemiLepEventHypothesis.jes  = float(lJesFactor)
process.hitFitTtSemiLepEventHypothesis.jesB = float(bJesFactor)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTop.root')
)

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::REDIGI311X", HLTPaths = ["HLT_Mu15_v*"], throw=True)

process.leadingJetSelection.src = 'tightLeadingPFJets'
process.bottomJetSelection.src  = 'tightBottomPFJets'

## end path   
process.path = cms.Path(#process.patDefaultSequence *
                        process.hltFilter *
                        process.semiLeptonicSelection *
                        process.semiLeptonicEvents *
                        process.makeGenEvt *
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
prependPF2PATSequence(process, options = {'runOnOLDcfg': True, 'runOnMC': True, 'electronIDs': '', 'switchOffEmbedding': False, 'skipIfNoPFMuon': True})

## adaptions (re-aranging of modules) to speed up processing
pathnames = process.paths_().keys()
for pathname in pathnames:
    ## move the ttGenEvent filter to the beginning of the sequence
    getattr(process, pathname).remove(process.genEvt)
    getattr(process, pathname).insert(0,process.genEvt)
    getattr(process, pathname).remove(process.decaySubset)
    getattr(process, pathname).insert(0,process.decaySubset)
    getattr(process, pathname).remove(process.initSubset)
    getattr(process, pathname).insert(0,process.initSubset)
    ## move the trigger to the beginning of the sequence
    getattr(process, pathname).remove(process.hltFilter)
    getattr(process, pathname).insert(0,process.hltFilter)


