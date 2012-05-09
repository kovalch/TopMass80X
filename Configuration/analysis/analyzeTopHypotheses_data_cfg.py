import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

# create lepton channel label 
options.register('lepton', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen decay channel")

# define the syntax for parsing
# you need to enter in the cfg file:
# search for arguments entered after cmsRun
if( hasattr(sys, "argv") ):
    # split arguments by comma - seperating different variables
    for args in sys.argv :
        arg = args.split(',')
        # split further by = to separate variable name and value
        for val in arg:
            val = val.split('=')
            # set variable var to value val (expected crab syntax: var=val)
            if(len(val)==2):
                setattr(options,val[0], val[1])

## choose the semileptonic decay channel (electron or muon)
#decayChannel=options.lepton
if(options.lepton=='unset'): 
    if(not globals().has_key('decayChannel')):
        decayChannel = 'muon' # 'electron'
else:
    decayChannel=options.lepton
print "used lepton decay channel: "+decayChannel

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

if (decayChannel=='electron'):
  readFiles.extend( [
         '/store/data/Run2011A/ElectronHad/AOD/PromptReco-v4/000/165/993/62F23D23-288B-E011-8F19-003048F1BF68.root'
  ] )
else:
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

## JET selection
process.tightBottomSSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") > 1.74'
                                           )
process.tightBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > 0.679'
                                           )

process.leadingJetSelection.src = 'tightLeadingPFJets'
process.bottomJetSelection.src  = 'tightBottomPFJets'

## b-tag selection
process.tightBottomPFJets.cut = 'bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679';
process.bottomJetSelection.minNumber = 2;

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "leps", "tightMuons")
if (decayChannel=='electron'):
  useElectronsForAllTtSemiLepHypotheses(process, 'goodElectronsEJ')
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "patMETsPF")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "combinedSecondaryVertexBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 0.679
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 0.679
process.hitFitTtSemiLepEventHypothesis.useBTagging       = True

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
                        process.tightBottomSSVPFJets *
                        process.tightBottomCSVPFJets *
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


## switch to from muon to electron collections
if (decayChannel=="electron"):
    process.TFileService.fileName = "analyzeTop.root"
    '/store/data/Run2011A/ElectronHad/AOD/PromptReco-v4/000/165/993/62F23D23-288B-E011-8F19-003048F1BF68.root'
    # adpat trigger
    process.hltFilter.HLTPaths = ["HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v*", "HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v3", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v*"]
    ## lepton-jet veto
    from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
    from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
    process.noOverlapJetsPFelec = cleanPatJets.clone(
        src = cms.InputTag("selectedPatJets"),
        preselection = cms.string(''),
        checkOverlaps = cms.PSet(
          electrons = cms.PSet(
            src       = cms.InputTag("tightElectronsEJ"),
            algorithm = cms.string("byDeltaR"),
            preselection        = cms.string(''),
            deltaR              = cms.double(0.3),
            checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            pairCut             = cms.string(""),
            requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
            )
          ),
        finalCut = cms.string(''),
        )
    process.goodJetsPF20.src  ='noOverlapJetsPFelec'
    process.centralJetsPF.src ='noOverlapJetsPFelec'
    process.reliableJetsPF.src='noOverlapJetsPFelec'
    process.noEtaJetsPF.src   ='noOverlapJetsPFelec'
    process.noPtJetsPF.src    ='noOverlapJetsPFelec'
    process.noConstJetsPF.src ='noOverlapJetsPFelec'
    process.noCEFJetsPF.src   ='noOverlapJetsPFelec'
    process.noNHFJetsPF.src   ='noOverlapJetsPFelec'
    process.noNEFJetsPF .src  ='noOverlapJetsPFelec'
    process.noCHFJetsPF.src   ='noOverlapJetsPFelec'
    process.noNCHJetsPF.src   ='noOverlapJetsPFelec'
    
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    pathnames = process.paths_().keys()
    for pathname in pathnames:
      # replace jet lepton veto
      getattr(process, pathname).replace(process.noOverlapJetsPF, process.noOverlapJetsPFelec)
      # replace muon selection
      getattr(process, pathname).replace(process.muonSelection, process.electronSelection)
      getattr(process, pathname).remove(process.secondMuonVeto)
      getattr(process, pathname).remove(process.electronVeto)
      # replace muon by electron in (remaining) kinfit analyzers
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'tightMuons', 'goodElectronsEJ')
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'noOverlapJetsPF', 'noOverlapJetsPFelec')

  
from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
prependPF2PATSequence(process, options = {'runOnOLDcfg': True,
                                          'runOnMC': False,
                                          'runOnAOD': True,
                                          'electronIDs': ['CiC','classical'],
                                          'switchOffEmbedding': False,
                                          'pfIsoConeMuon': 0.4,
                                          'pfIsoConeElec': 0.4,
                                          #'skipIfNoPFMuon': True,
                                          'METCorrectionLevel': 2,
                                          })

## adaptions (re-aranging of modules) to speed up processing
pathnames = process.paths_().keys()
for pathname in pathnames:
    ## move the trigger to the beginning of the sequence
    getattr(process, pathname).remove(process.hltFilter)
    getattr(process, pathname).insert(0,process.hltFilter)


