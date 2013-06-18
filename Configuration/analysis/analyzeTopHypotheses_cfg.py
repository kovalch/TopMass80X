import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

 # for Summer11/Fall11 MC one can choose: ttbar, wjets, zjets, singleAntiTopS, singleTopT, singleAntiTopT, singleTopTw, singleAntiTopTw, singleTopS WW, WZ, ZZ, qcd (for muon channel); qcdEM1, qcdEM2, qcdEM3, qcdBCE1, qcdBCE2, qcdBCE3 (for electron channel), zprime_m500gev_w5000mev, zprime_m750gev_w7500mev
options.register('sample', 'none',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen sample")
options.register('mcversion', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "mcversion of chosen sample")
# create lepton channel label 
options.register('lepton', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen decay channel")
options.register('metcl', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "MET correction level")
options.register('pdf', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Calculate PDF weights")


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

## version of the MC (important for trigger)
## This parameter is also used to select the correct procedure for PU event reweighting
if(options.mcversion=='unset'): 
  MCversion = 'Summer11'
else:
  MCversion = options.mcversion
print "used mcversion: "+MCversion

## choose the semileptonic decay channel (electron or muon)
#decayChannel=options.lepton
if(options.lepton=='unset'): 
    if(not globals().has_key('decayChannel')):
        decayChannel = 'muon' # 'electron'
else:
    decayChannel=options.lepton
print "used lepton decay channel: "+decayChannel

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

mesShift = '@mesShift@'
if mesShift.startswith('@'):
  mesShift = '0.0'
if (mesShift == 'down'):
  mesShift = '-1.0'
if (mesShift == 'up'):
  mesShift = '1.0'

print mesShift

eesShift = '@eesShift@'
if eesShift.startswith('@'):
  eesShift = '0.0'
if (eesShift == 'down'):
  eesShift = '-1.0'
if (eesShift == "up"):
  eesShift = '1.0'

uncFactor = '@uncFactor@'
if uncFactor.startswith('@'):
  uncFactor = '1.0'

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
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,
                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
                              inputCommands=cms.untracked.vstring(
                                      'keep *',
                                      'drop LHERunInfoProduct_*_*_*'
                              )
                              #,skipEvents=cms.untracked.uint32(4000)
)
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  readFiles.extend( [
         '/store/user/eschliec/TTJets_TuneD6T_7TeV-madgraph-tauola/PATWithPF_v4/e59efddd8a1547799dca5b47d5556447/patTuple_9_1_Nb0.root'
  ] )
if (MCversion == "Fall11"):
  readFiles.extend( [
         '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/00CAD0AC-17FA-E011-A4F4-00304867905A.root'
  ] )
if (MCversion == "Summer12"):
  readFiles.extend( [
         '/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/FC5CECAE-8B14-E211-8578-0025B3E0652A.root'
  ] )
else:
  readFiles.extend( [
  #       'file:/scratch/hh/current/cms/user/stadie/2011/production/TT_TuneZ2_7TeV_pythia6_FASTSIM/TT_TuneZ2_7TeV_pythia6_FASTSIM_0.root',
  #       '/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/00454769-577B-E011-ACCD-001E0B49808A.root',
         #'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FEEE3638-F297-E011-AAF8-00304867BEC0.root',
  #       '/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/02719D6B-1398-E011-AA71-001A92971B94.root',
  #       '/store/mc/Summer11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0004EB5E-64AC-E011-B046-003048678FDE.root',
  #       '/store/mc/Fall11/TT_TuneZ2_7TeV-mcatnlo/AODSIM/PU_S6_START42_V14B-v1/0001/FED956DF-7F2A-E111-8124-002618943958.root'
          '/store/user/mgosseli/mc/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1_139.root',
          '/store/user/mgosseli/mc/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1_683.root'
  ] )

secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

## configure geometry & conditions
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.GlobalTag.globaltag = cms.string('START41_V0::All')
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.GlobalTag.globaltag = cms.string('START42_V17::All')
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
  process.GlobalTag.globaltag = cms.string('START53_V7A::All')

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

## electron shift
process.load("TopAnalysis.TopUtils.ElectronEnergyScale_cfi")
from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *

process.scaledElectronEnergy.src      = "selectedPatElectrons"
process.scaledElectronEnergy.mets     = "scaledJetEnergy:patMETs"
process.scaledElectronEnergy.shiftBy  = float(eesShift)
process.vertexSelectedElectrons.src   = "scaledElectronEnergy:selectedPatElectrons"

## muon shift
process.load("TopAnalysis.TopUtils.MuonEnergyScale_cfi")
from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *

process.scaledMuonEnergy.src      = "selectedPatMuons"
process.scaledMuonEnergy.mets     = "scaledElectronEnergy:METs"
process.scaledMuonEnergy.shiftBy  = float(mesShift)
process.vertexSelectedMuons.src   = "scaledMuonEnergy:selectedPatMuons"

## unclustered energy scale
process.load("TopAnalysis.TopUtils.UnclusteredMETScale_cfi")
from TopAnalysis.TopUtils.UnclusteredMETScale_cfi import *

process.scaledMET.inputJets       = "scaledJetEnergy:selectedPatJets"
process.scaledMET.inputMETs       = "scaledMuonEnergy:METs"
process.scaledMET.inputElectrons  = "scaledElectronEnergy:selectedPatElectrons"
process.scaledMET.inputMuons      = "scaledMuonEnergy:selectedPatMuons"
process.scaledMET.scaleFactor     = float(uncFactor)

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
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu24_eta2p1_v*"], throw=True)
if(MCversion=="Summer11"):
    process.hltFilter.HLTPaths=["HLT_IsoMu24_v*"]

## JET selection
process.tightBottomPFJets.cut = 'bDiscriminator("combinedSecondaryVertexBJetTags") > 0.5';

process.tightBottomSSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") > 1.74'
                                           )
process.tightBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > 0.679'
                                           )
process.looseBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > 0.5'
                                           )

process.leadingJetSelection.src = 'tightLeadingPFJets'

process.bottomJetSelection.src  = 'looseBottomCSVPFJets'
process.bottomJetSelection.minNumber = 2;

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "leps", "tightMuons")
if (decayChannel=='electron'):
  useElectronsForAllTtSemiLepHypotheses(process, 'goodElectronsEJ')
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "scaledMET:scaledMETs")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
process.ttSemiLepJetPartonMatch.maxNJets  = -1

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "combinedSecondaryVertexBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 0.5
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 0.85
process.hitFitTtSemiLepEventHypothesis.useBTagging       = True

addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )
#removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
process.load("TopMass.TopEventTree.EventHypothesisAnalyzer_cff")

## PDF weights
#process.analyzeHitFit.savePDFWeights = options.pdf

## MC weights
process.load("TopAnalysis.TopUtils.EventWeightMC_cfi")

## ============================
##  MC PU reweighting
## ============================

process.load("TopAnalysis.TopUtils.EventWeightPU_cfi")

## Apply common setting before module is cloned for systematic studies

process.eventWeightPU.MCSampleTag = MCversion

if (MCversion == "Fall11"):
    process.eventWeightPU.MCSampleHistoName        = "histo_Fall11_true"
    process.eventWeightPU.DataHistoName            = "histoData_true"
elif (MCversion == "Summer11"):    
    process.eventWeightPU.MCSampleHistoName        = "histoSummer11_flat_true"
    process.eventWeightPU.DataHistoName            = "histoData_true_fineBinning"

process.eventWeightPUsysNo   = process.eventWeightPU.clone()
process.eventWeightPUsysUp   = process.eventWeightPU.clone()
process.eventWeightPUsysDown = process.eventWeightPU.clone()

#### Parameters 'CreateWeight3DHisto' and 'Weight3DHistoFile' required for cff-file, but actually not used for Fall11 samples
    
#### Configuration for Nominal PU Weights

process.eventWeightPUsysNo.WeightName          = "eventWeightPU"
process.eventWeightPUsysNo.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysNo_68000_2011Full.root"
process.eventWeightPUsysNo.CreateWeight3DHisto = False
#process.eventWeightPUsysNo.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3D.root"

#### Configuration for PU Up Variations

process.eventWeightPUsysUp.WeightName          = "eventWeightPUUp"
process.eventWeightPUsysUp.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_71400_2011Full.root"
process.eventWeightPUsysUp.CreateWeight3DHisto = False
#process.eventWeightPUsysUp.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3DUp.root"

#### Configuration for PU Down Variations

process.eventWeightPUsysDown.WeightName          = "eventWeightPUDown"
process.eventWeightPUsysDown.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysDown_64600_2011Full.root"
process.eventWeightPUsysDown.CreateWeight3DHisto = False
#process.eventWeightPUsysDown.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3DDown.root"

#### event weight sequence

process.makeEventWeightsPU = cms.Sequence(process.eventWeightPUsysNo   *
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
process.bTagSFEventWeight.bTagAlgo = "CSVM"
process.bTagSFEventWeight.version  = "11-004"
process.bTagSFEventWeight.sysVar   = "" # bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown possible;
process.bTagSFEventWeight.filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiencyCSVM.root"
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
process.effSFMuonEventWeight.filename              = "TopAnalysis/Configuration/data/MuonEffSF2011.root"
process.effSFMuonEventWeight.verbose               = 0
process.effSFMuonEventWeight.additionalFactor      = 1. ## lepton selection and trigger eff. SF both included in loaded histo
process.effSFMuonEventWeight.additionalFactorErr   = 0.01 ## 1% sys error to account for non-flatness
process.effSFMuonEventWeight.meanTriggerEffSF      = 0.9824
process.effSFMuonEventWeight.shapeDistortionFactor = -1

process.load("TopAnalysis.TopUtils.EffSFElectronEventWeight_cfi")
process.effSFElectronEventWeight.electrons=cms.InputTag("goodElectronsEJ")
process.effSFElectronEventWeight.jets=cms.InputTag("tightLeadingPFJets")
process.effSFElectronEventWeight.sysVar   = cms.string("")
process.effSFElectronEventWeight.verbose=cms.int32(0)
process.effSFElectronEventWeight.filenameJetLeg="TopAnalysis/Configuration/data/JetLegTriggerEfficiencyIsoLepTriJetJetMult4.root"
process.effSFElectronEventWeight.additionalFactor=1. ## lepton selection eff. SF
process.effSFElectronEventWeight.additionalFactorErr=0.02 ## 2% sys error to account for selection difference Z - ttbar
process.effSFElectronEventWeight.meanTriggerEffSF=0.968
process.effSFElectronEventWeight.meanTriggerEffSFErr=0.004
process.effSFElectronEventWeight.shapeDistortionErr=0.02
process.effSFElectronEventWeight.jetTriggerEffsSFNormSysErr =0.01
process.effSFElectronEventWeight.jetTriggerEffsSFShapeSysErr=0.005


# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTop.root')
)

# register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = "Tree for UHH top-quark analysis\nParticles are in order {TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino}"

#process.MessageLogger = cms.Service("MessageLogger")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

## end path   
process.path = cms.Path(#process.patDefaultSequence *
                        process.hltFilter *
                        process.scaledJetEnergy *
                        process.scaledElectronEnergy *
                        #process.content *
                        process.scaledMuonEnergy *
                        #process.content *
                        process.scaledMET *
                        process.semiLeptonicSelection *
                        process.tightBottomSSVPFJets *
                        process.tightBottomCSVPFJets *
                        process.looseBottomCSVPFJets *
                        process.semiLeptonicEvents *
                        #process.eventWeightMC *
                        #process.makeEventWeightsPU *
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
#process.path.remove(process.unconvTightElectronsEJ)
#process.path.remove(process.goodElectronsEJ)
#process.path.remove(process.looseElectronsEJ)
#process.path.remove(process.tightElectronsEJ)



## switch to from muon to electron collections
if (decayChannel=="electron"):
    process.TFileService.fileName = "analyzeTop.root"
    # adpat trigger
    if (MCversion == "Fall11"):
        process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
    else:
        process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v*"]
    ## lepton-jet veto
    from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
    from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
    process.noOverlapJetsPFelec = cleanPatJets.clone(
        src = cms.InputTag("scaledJetEnergy:selectedPatJets"),
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
      ## replace effSF
      getattr(process, pathname).replace(process.effSFMuonEventWeight, process.effSFElectronEventWeight)
      # replace muon by electron in (remaining) kinfit analyzers
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'tightMuons', 'goodElectronsEJ')
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'noOverlapJetsPF', 'noOverlapJetsPFelec')
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'effSFMuonEventWeight', 'effSFElectronEventWeight')



from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
prependPF2PATSequence(process, options = {'runOnOLDcfg': True,
                                          'runOnMC': True,
                                          'runOnAOD': True,
                                          'electronIDs': ['CiC','classical','MVA'],
                                          'switchOffEmbedding': False,
                                          'pfIsoConeMuon': 0.4,
                                          'pfIsoConeElec': 0.4,
                                          #'skipIfNoPFMuon': True,
                                          'METCorrectionLevel': options.metcl,
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


