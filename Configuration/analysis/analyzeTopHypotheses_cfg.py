import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

# for summer11 MC one can choose: ttbar, wjets, zjets, singleAntiTopS, singleTopT, singleAntiTopT, singleTopTw, singleAntiTopTw, WW, WZ, qcd (for muon channel);
# still missing: ZZ, singleTopS
options.register('sample', 'none',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen sample")

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

# print chosen sample (e.g. ttbar)
# value is known from external parsing
# if set, switches runOnAOD in PF2PAT to true
print "Chosen sample to run over: ", options.sample

## top mass measurement
process = cms.Process("topMass")

## grid control parameters
lJesFactor = '@lJesFactor@'
if lJesFactor.startswith('@'):
  lJesFactor = '1.0'

bJesFactor = '@bJesFactor@'
if bJesFactor.startswith('@'):
  bJesFactor = '1.0'

resolution = '@resolution@'
if resolution.startswith('@'):
  resolution = 'nominal'

scaleType = '@scaleType@'
if scaleType.startswith('@'):
  scaleType = 'abs'

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
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  readFiles.extend( [
         '/store/user/eschliec/TTJets_TuneD6T_7TeV-madgraph-tauola/PATWithPF_v4/e59efddd8a1547799dca5b47d5556447/patTuple_9_1_Nb0.root'
  ] )
if(options.sample=="ttbarFall11"):
  readFiles.extend( [
         '/store/user/eschliec/TTJets_TuneD6T_7TeV-madgraph-tauola/PATWithPF_v4/e59efddd8a1547799dca5b47d5556447/patTuple_9_1_Nb0.root'
  ] )
else:
  readFiles.extend( [
  #       'rfio:/scratch/hh/current/cms/user/stadie/AODIntegrationTestWithHLT.root',
  #       '/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/00454769-577B-E011-ACCD-001E0B49808A.root',
         '/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FEEE3638-F297-E011-AAF8-00304867BEC0.root',
  #       '/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/02719D6B-1398-E011-AA71-001A92971B94.root',
  #       '/store/mc/Summer11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0004EB5E-64AC-E011-B046-003048678FDE.root',
  ] )

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
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.GlobalTag.globaltag = cms.string('START41_V0::All')
else:
  process.GlobalTag.globaltag = cms.string('START42_V17::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")

## redefine veto jets to be sure it is also replaced when running on PF
from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJets
process.vetoJets.src="goodJetsPF30"
process.vetoJets.cut=''

## configure JetEnergyScale tool
process.load("TopAnalysis.TopUtils.JetEnergyScale_cff")
from TopAnalysis.TopUtils.JetEnergyScale_cff import *

scaledJetEnergy.scaleType    = scaleType
scaledJetEnergy.inputJets    = "selectedPatJetsAK5PF"
scaledJetEnergy.inputMETs    = "patMETsPF"
scaledJetEnergy.scaleFactor  = float(lJesFactor)
scaledJetEnergy.scaleFactorB = float(bJesFactor)
if (resolution=='down'):
  scaledJetEnergy.resolutionFactors   = [0.994, 1.126, 1.006, 0.961]
if (resolution=='nominal'):
  scaledJetEnergy.resolutionFactors   = [1.066, 1.191, 1.096, 1.166]
if (resolution=='up'):
  scaledJetEnergy.resolutionFactors   = [1.140, 1.258, 1.190, 1.370]
scaledJetEnergy.resolutionEtaRanges   = [0, 1.1, 1.1, 1.7, 1.7, 2.3, 2.3, -1]


process.noOverlapJetsPF.src = "scaledJetEnergy:selectedPatJets"

## sequences for ttGenEvent and TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

## TRIGGER
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::REDIGI311X", HLTPaths = ["HLT_Mu15_v*"], throw=True)
else:
  process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu17_v*"], throw=True)

## JET selection
#process.tightBottomPFJets.cut = 'bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74';

process.tightBottomSSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") > 1.74'
                                           )
process.tightBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > 0.679'
                                           )

process.leadingJetSelection.src = 'tightLeadingPFJets'

process.bottomJetSelection.src  = 'tightBottomPFJets'
process.bottomJetSelection.minNumber = 0;

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "leps", "tightMuons")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "scaledJetEnergy:patMETs")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "simpleSecondaryVertexHighEffBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 1.74
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 1.74
process.hitFitTtSemiLepEventHypothesis.useBTagging       = False

addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )
#removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
process.load("TopMass.Analyzer.EventHypothesisAnalyzer_cff")

## ============================
##  MC PU reweighting
## ============================

process.load("TopAnalysis.TopUtils.EventWeightPU_cfi")

process.eventWeightPU        = process.eventWeightPU.clone()
process.eventWeightPUsysUp   = process.eventWeightPU.clone()
process.eventWeightPUsysDown = process.eventWeightPU.clone()

#### Configuration for Nominal PU Weights

process.eventWeightPU.WeightName          = "eventWeightPU"
process.eventWeightPU.Weight3DName        = "eventWeightPU3D"
process.eventWeightPU.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_2011Full.root"
process.eventWeightPU.Data3DFile          = "TopAnalysis/TopUtils/data/Data_PUDist_2011Full.root"

process.eventWeightPU.CreateWeight3DHisto = False
process.eventWeightPU.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3D.root"

#### Configuration for PU Up Variations

process.eventWeightPUsysUp.WeightName          = "eventWeightPUUp"
process.eventWeightPUsysUp.Weight3DName        = "eventWeightPU3DUp"
process.eventWeightPUsysUp.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_2011Full.root"
process.eventWeightPUsysUp.Data3DFile          = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_2011Full.root"

process.eventWeightPUsysUp.CreateWeight3DHisto = False
process.eventWeightPUsysUp.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3DUp.root"

#### Configuration for PU Down Variations

process.eventWeightPUsysDown.WeightName          = "eventWeightPUDown"
process.eventWeightPUsysDown.Weight3DName        = "eventWeightPU3DDown"
process.eventWeightPUsysDown.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysDown_2011Full.root"
process.eventWeightPUsysDown.Data3DFile          = "TopAnalysis/TopUtils/data/Data_PUDist_sysDown_2011Full.root"

process.eventWeightPUsysDown.CreateWeight3DHisto = False
process.eventWeightPUsysDown.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3DDown.root"

process.makeEventWeightsPU = cms.Sequence(process.eventWeightPU        *
                                          process.eventWeightPUsysUp   *
                                          process.eventWeightPUsysDown  )

## ---
##    MC B-tag reweighting
## ---
## load BTV database
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

process.load("TopAnalysis.TopUtils.BTagSFEventWeight_cfi")
process.bTagSFEventWeight.jets     = "tightLeadingPFJets"
process.bTagSFEventWeight.bTagAlgo = "SSVHEM"
process.bTagSFEventWeight.sysVar   = "" # bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown possible;
process.bTagSFEventWeight.filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiency.root"
process.bTagSFEventWeight.verbose  = 0

process.bTagSFEventWeightBTagSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
process.bTagSFEventWeightBTagSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
process.bTagSFEventWeightMisTagSFUp   = process.bTagSFEventWeight.clone(sysVar = "misTagSFUp")
process.bTagSFEventWeightMisTagSFDown = process.bTagSFEventWeight.clone(sysVar = "misTagSFDown")


## ---
##    MC eff SF reweighting
## ---
## scale factor for trigger and lepton selection efficiency
process.load("TopAnalysis.TopUtils.EffSFMuonEventWeight_cfi")
process.effSFMuonEventWeight.particles             = "tightMuons"
process.effSFMuonEventWeight.sysVar                = ""
process.effSFMuonEventWeight.filename              = "TopAnalysis/Configuration/data/efficiencyIsoMu17Combined_tapTrigger_SF_Eta.root"
process.effSFMuonEventWeight.verbose               = 0
process.effSFMuonEventWeight.additionalFactor      = 0.9990 ## lepton selection eff. SF
process.effSFMuonEventWeight.additionalFactorErr   = 0.03 ## 3% sys error to account for selection difference Z - ttbar
process.effSFMuonEventWeight.meanTriggerEffSF      = 0.9905
process.effSFMuonEventWeight.shapeDistortionFactor = 0.5



# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTop.root')
)

## end path   
process.path = cms.Path(#process.patDefaultSequence *
                        process.hltFilter *
                        process.scaledJetEnergy *
                        process.semiLeptonicSelection *
                        process.tightBottomSSVPFJets *
                        process.tightBottomCSVPFJets *
                        process.semiLeptonicEvents *
                        process.makeEventWeightsPU *
                        process.bTagSFEventWeight *
                        process.bTagSFEventWeightBTagSFUp *
                        process.bTagSFEventWeightBTagSFDown *
                        process.bTagSFEventWeightMisTagSFUp *
                        process.bTagSFEventWeightMisTagSFDown *
                        process.effSFMuonEventWeight *
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
prependPF2PATSequence(process, options = {'runOnOLDcfg': True,
                                          'runOnMC': True,
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


