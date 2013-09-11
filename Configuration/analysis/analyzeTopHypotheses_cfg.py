import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('mcversion', 'Summer12', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC campaign or data")
options.register('mcWeight', 1.0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "MC sample event weight")
options.register('lepton', 'muon', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Lepton+jets channel")
options.register('metcl', 1, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "MET correction level")

options.register('scaleType', 'abs', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "JES scale type")
options.register('jessource', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Uncertainty source for JES variation for source:up/down")
options.register('lJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "JES")
options.register('bJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "bJES")
options.register('resolution', 'nominal', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "JER")
options.register('mesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Muon energy scale shift in sigma")
options.register('eesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Electron energy scale shift in sigma")
options.register('uncFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Unclustered energy factor")

options.register('csvm', 0.679, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "CSVM working point")
options.register('nbjets', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Minimum number of bjets")

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

print options

if (options.mcversion == "data"): data = True
else:                             data = False

## top mass measurement
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
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,
                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
                              inputCommands=cms.untracked.vstring(
                                      'keep *',
                                      'drop LHERunInfoProduct_*_*_*'
                              )
)
if (options.mcversion == "Fall11"):
  readFiles.extend( [
         '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/00CAD0AC-17FA-E011-A4F4-00304867905A.root'
  ] )
elif (options.mcversion == "Summer12"):
  readFiles.extend( [
         '/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/FC5CECAE-8B14-E211-8578-0025B3E0652A.root'
  ] )
elif (data and options.lepton == "muon"):
  readFiles.extend( [
         '/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/30000/FEDCB8E2-5270-E211-8FD6-00266CFFBC38.root'
  ] )
elif (data and options.lepton == "electron"):
  readFiles.extend( [
         '/store/data/Run2012A/SingleElectron/AOD/22Jan2013-v1/30002/72352080-9372-E211-B431-00266CFFA1AC.root'
  ] )
  process.source.skipEvents=cms.untracked.uint32(1800)
else:
  readFiles.extend( [
          '/store/user/mgosseli/mc/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1_139.root',
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
if data:
  process.GlobalTag.globaltag = cms.string('FT_53_V21_AN5::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")

## redefine veto jets to be sure it is also replaced when running on PF
from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJets
process.vetoJets.src = "goodJetsPF30"
process.vetoJets.cut = ''

if not data:
    ## configure JetEnergyScale tool
    process.load("TopAnalysis.TopUtils.JetEnergyScale_cff")
    from TopAnalysis.TopUtils.JetEnergyScale_cff import *

    scaledJetEnergy.scaleType    = options.scaleType
    scaledJetEnergy.JECUncSrcFile= "TopAnalysis/TopUtils/data/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt"
    scaledJetEnergy.sourceName   = options.jessource
    scaledJetEnergy.scaleFactor  = options.lJesFactor
    scaledJetEnergy.scaleFactorB = options.bJesFactor
    if (options.resolution=='down'):
      scaledJetEnergy.resolutionFactors   = [0.990 , 1.001 , 1.032 , 1.042 , 1.089]
    if (options.resolution=='nominal'):
      scaledJetEnergy.resolutionFactors   = [1.052 , 1.057 , 1.096 , 1.134 , 1.288]
    if (options.resolution=='up'):
      scaledJetEnergy.resolutionFactors   = [1.115 , 1.114 , 1.161 , 1.228 , 1.488]
    scaledJetEnergy.resolutionEtaRanges   = [0.0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,-1.]

    scaledJetEnergy.inputJets    = "selectedPatJets"
    scaledJetEnergy.inputMETs    = "patMETs"

    process.goodJetsPF.src = "scaledJetEnergy:selectedPatJets"

    ## electron shift
    process.load("TopAnalysis.TopUtils.ElectronEnergyScale_cfi")
    from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *

    process.scaledElectronEnergy.src      = "selectedPatElectrons"
    process.scaledElectronEnergy.mets     = "scaledJetEnergy:patMETs"
    process.scaledElectronEnergy.shiftBy  = options.eesShift
    process.vertexSelectedElectrons.src   = "scaledElectronEnergy:selectedPatElectrons"

    ## muon shift
    process.load("TopAnalysis.TopUtils.MuonEnergyScale_cfi")
    from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *

    process.scaledMuonEnergy.src      = "selectedPatMuons"
    process.scaledMuonEnergy.mets     = "scaledElectronEnergy:METs"
    process.scaledMuonEnergy.shiftBy  = options.mesShift
    process.vertexSelectedMuons.src   = "scaledMuonEnergy:selectedPatMuons"

    ## unclustered energy scale
    process.load("TopAnalysis.TopUtils.UnclusteredMETScale_cfi")
    from TopAnalysis.TopUtils.UnclusteredMETScale_cfi import *

    process.scaledMET.inputJets       = "scaledJetEnergy:selectedPatJets"
    process.scaledMET.inputMETs       = "scaledMuonEnergy:METs"
    process.scaledMET.inputElectrons  = "scaledElectronEnergy:selectedPatElectrons"
    process.scaledMET.inputMuons      = "scaledMuonEnergy:selectedPatMuons"
    process.scaledMET.scaleFactor     = options.uncFactor

    ## sequence for ttGenEvent
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## sequence for TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

## TRIGGER
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::REDIGI311X", HLTPaths = ["HLT_Mu15_v*"], throw=True)
else:
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu24_eta2p1_v*"], throw=True)
if(options.mcversion=="Summer11"):
    process.hltFilter.HLTPaths=["HLT_IsoMu24_v*"]
if(options.mcversion=="Summer12"):
    process.hltFilter.HLTPaths=["HLT_IsoMu24_eta2p1_v*"]
if data:
    process.hltFilter.HLTPaths=["HLT_IsoMu24_eta2p1_v*"]


## JET selection
process.tightBottomSSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") > 1.74'
                                           )
process.tightBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
                                           cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > ' + str(options.csvm)
                                           )

process.leadingJetSelection.src = 'tightLeadingPFJets'

process.bottomJetSelection.src  = 'tightBottomCSVPFJets'
process.bottomJetSelection.minNumber = options.nbjets;

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "leps", "tightMuons")
if (options.lepton=='electron'):
    useElectronsForAllTtSemiLepHypotheses(process, 'goodElectronsEJ')
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "scaledMET:scaledMETs")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
if data:
    setForAllTtSemiLepHypotheses(process, "mets", "patMETsPF")
    setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
process.ttSemiLepJetPartonMatch.maxNJets  = -1

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "combinedSecondaryVertexBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = options.csvm
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = options.csvm
process.hitFitTtSemiLepEventHypothesis.useBTagging       = True

addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )
if data: removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import analyzeHypothesis
process.analyzeHitFit = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key")
from TopMass.TopEventTree.JetEventAnalyzer_cfi import analyzeJets
process.analyzeJets = analyzeJets.clone()
from TopMass.TopEventTree.WeightEventAnalyzer_cfi import analyzeWeights
process.analyzeWeights = analyzeWeights.clone(
                                              mcWeight        = options.mcWeight,
                                              puWeightSrc     = cms.InputTag("eventWeightPUsysNo"  , "eventWeightPU"),
                                              puWeightUpSrc   = cms.InputTag("eventWeightPUsysUp"  , "eventWeightPU"),
                                              puWeightDownSrc = cms.InputTag("eventWeightPUsysDown", "eventWeightPU"),
                                              savePDFWeights = True
                                             )

if not data:
    ## MC weights
    process.load("TopAnalysis.TopUtils.EventWeightMC_cfi")

    ## ============================
    ##  MC PU reweighting
    ## ============================

    process.load("TopAnalysis.TopUtils.EventWeightPU_cfi")

    ## Apply common setting before module is cloned for systematic studies

    process.eventWeightPU.MCSampleTag = options.mcversion

    if (options.mcversion == "Fall11"):
        process.eventWeightPU.MCSampleHistoName        = "histo_Fall11_true"
        process.eventWeightPU.DataHistoName            = "histoData_true"
    elif (options.mcversion == "Summer11"):    
        process.eventWeightPU.MCSampleHistoName        = "histoSummer11_flat_true"
        process.eventWeightPU.DataHistoName            = "histoData_true_fineBinning"
    elif (options.mcversion == "Summer12"):
        process.eventWeightPU.MCSampleHistoName        = "puhisto"
        process.eventWeightPU.DataHistoName            = "pileup"

        process.eventWeightPU.MCSampleFile             = "TopAnalysis/TopUtils/data/MC_PUDist_Summer12_S10.root"

    process.eventWeightPUsysNo   = process.eventWeightPU.clone()
    process.eventWeightPUsysUp   = process.eventWeightPU.clone()
    process.eventWeightPUsysDown = process.eventWeightPU.clone()

    #### Parameters 'CreateWeight3DHisto' and 'Weight3DHistoFile' required for cff-file, but actually not used for Fall11 samples
        
    #### Configuration for Nominal PU Weights

    process.eventWeightPUsysNo.WeightName          = "eventWeightPU"
    process.eventWeightPUsysNo.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysNo_69400_2012ABCD22JanReReco_190456-208686_8TeV.root"
    process.eventWeightPUsysNo.CreateWeight3DHisto = False
    #process.eventWeightPUsysNo.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3D.root"

    #### Configuration for PU Up Variations

    process.eventWeightPUsysUp.WeightName          = "eventWeightPUUp"
    process.eventWeightPUsysUp.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_73564_2012ABCD22JanReReco_190456-208686_8TeV.root"
    process.eventWeightPUsysUp.CreateWeight3DHisto = False
    #process.eventWeightPUsysUp.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3DUp.root"

    #### Configuration for PU Down Variations

    process.eventWeightPUsysDown.WeightName          = "eventWeightPUDown"
    process.eventWeightPUsysDown.DataFile            = "TopAnalysis/TopUtils/data/Data_PUDist_sysDn_65236_2012ABCD22JanReReco_190456-208686_8TeV.root"
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
    process.bTagSFEventWeight.version  = "2012"
    process.bTagSFEventWeight.sysVar   = "" # bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown possible;
    process.bTagSFEventWeight.filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiency2012.root" 
    process.bTagSFEventWeight.verbose  = 0

    process.bTagSFEventWeightBTagSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
    process.bTagSFEventWeightBTagSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
    process.bTagSFEventWeightMisTagSFUp   = process.bTagSFEventWeight.clone(sysVar = "misTagSFUp")
    process.bTagSFEventWeightMisTagSFDown = process.bTagSFEventWeight.clone(sysVar = "misTagSFDown")


    ## ---
    ##    MC eff SF reweighting
    ## ---
    ## scale factor for trigger and lepton selection efficiency
    process.load("TopAnalysis.TopUtils.EffSFLepton2DEventWeight_cfi")
    process.effSFMuonEventWeight=process.effSFLepton2DEventWeight.clone()
    process.effSFMuonEventWeight.particles=cms.InputTag("tightMuons")
    process.effSFMuonEventWeight.sysVar   = cms.string("")
    process.effSFMuonEventWeight.filename= "TopAnalysis/Configuration/data/MuonEffSF2D2012.root"
    process.effSFMuonEventWeight.verbose=cms.int32(0)
    process.effSFMuonEventWeight.additionalFactor=1.0 ## lepton selection and trigger eff. SF both included in loaded histo
    process.effSFMuonEventWeight.additionalFactorErr=0.01 ## 1.0% sys error to account for selection difference Z - ttbar
    process.effSFMuonEventWeight.shapeDistortionFactor=-1

    process.effSFElectronEventWeight=process.effSFLepton2DEventWeight.clone()
    process.effSFElectronEventWeight.particles=cms.InputTag("goodElectronsEJ")
    process.effSFElectronEventWeight.jets=cms.InputTag("tightLeadingPFJets")
    process.effSFElectronEventWeight.sysVar   = cms.string("")
    process.effSFElectronEventWeight.filename= "TopAnalysis/Configuration/data/EleEffSF2D2012.root"
    process.effSFElectronEventWeight.verbose=cms.int32(0)
    process.effSFElectronEventWeight.additionalFactor=1.0 ## lepton selection and trigger eff. SF both included in loaded histo
    process.effSFElectronEventWeight.additionalFactorErr=0.01 ## 1.0% sys error to account for selection difference Z - ttbar
    process.effSFElectronEventWeight.shapeDistortionFactor=-1


# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTop.root')
)

# register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = ""

#process.MessageLogger = cms.Service("MessageLogger")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

## end path   
if data:
    process.path = cms.Path(
                            process.hltFilter *
                            process.semiLeptonicSelection *
                            process.tightBottomSSVPFJets *
                            process.tightBottomCSVPFJets *
                            process.semiLeptonicEvents *
                            process.makeTtSemiLepEvent *
                            process.analyzeHitFit *
                            process.analyzeJets *
                            process.analyzeWeights
                            )
else:
    process.path = cms.Path(
                            process.hltFilter *
                            process.scaledJetEnergy *
                            process.scaledElectronEnergy *
                            process.scaledMuonEnergy *
                            process.scaledMET *
                            process.semiLeptonicSelection *
                            process.tightBottomSSVPFJets *
                            process.tightBottomCSVPFJets *
                            process.semiLeptonicEvents *
                            process.eventWeightMC *
                            process.makeEventWeightsPU *
                            process.bTagSFEventWeight *
                            process.bTagSFEventWeightBTagSFUp *
                            process.bTagSFEventWeightBTagSFDown *
                            process.bTagSFEventWeightMisTagSFUp *
                            process.bTagSFEventWeightMisTagSFDown *
                            process.effSFMuonEventWeight *
                            process.makeGenEvt *
                            process.makeTtSemiLepEvent *
                            process.analyzeHitFit *
                            process.analyzeJets *
                            process.analyzeWeights
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
if (options.lepton=="electron"):
    process.analyzeHitFit.lepton = 11
    process.TFileService.fileName = "analyzeTop.root"
    # adpat trigger
    if (options.mcversion == "Fall11"):
        process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
    elif (options.mcversion == "Summer12"):
        process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
    elif data:
        process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
    ## lepton-jet veto
    from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
    from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
    process.noOverlapJetsPFelec = cleanPatJets.clone(
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
    if not data: process.noOverlapJetsPFelec.src = cms.InputTag("scaledJetEnergy:selectedPatJets")
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
      if not data: getattr(process, pathname).replace(process.effSFMuonEventWeight, process.effSFElectronEventWeight)
      # replace muon by electron in (remaining) kinfit analyzers
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'tightMuons', 'goodElectronsEJ')
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'effSFMuonEventWeight', 'effSFElectronEventWeight')



from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
PFoptions = {
        'runOnMC': not data,
        'runOnAOD': True,
        'switchOffEmbedding': False,
        'addResolutions': True,
        'resolutionsVersion': 'fall11',
        'runOnOLDcfg': True,
        'cutsMuon': 'pt > 10. & abs(eta) < 2.5',
        'cutsElec': 'et > 20. & abs(eta) < 2.5',
        'cutsJets': 'pt > 10 & abs(eta) < 5.0', 
        'electronIDs': ['CiC','classical','MVA'],
        'pfIsoConeMuon': 0.4,
        'pfIsoConeElec': 0.3,
        'pfIsoValMuon': 0.2,
        'pfIsoValElec': 0.15,
        'doDeltaBetaCorrMuon' : True,
        'doDeltaBetaCorrElec' : True,
        'skipIfNoPFMuon': False,
        'skipIfNoPFElec': False,
        'addNoCutPFMuon': False,
        'addNoCutPFElec': False,
        'noMuonTopProjection': False,
        'noElecTopProjection': False,
        'analyzersBeforeMuonIso':cms.Sequence(),
        'analyzersBeforeElecIso':cms.Sequence(),
        'excludeElectronsFromWsFromGenJets': True,
        'METCorrectionLevel': options.metcl,
        }
prependPF2PATSequence(process, options = PFoptions)

## adaptions (re-aranging of modules) to speed up processing
pathnames = process.paths_().keys()
for pathname in pathnames:
    if not data:
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

## adjust lepton pT cut
import re
ptval="33" # adjust cut value here! [GeV]
relevantLeptonCollections = [process.tightElectronsEJ, process.goldenMuons]
exp = re.compile('(?:t\s?>\s?30)') 
for lep in relevantLeptonCollections:
    if(exp.search(lep.cut.pythonValue())!=None):
        tmpExp=exp.sub("t > "+ptval, str(lep.cut.pythonValue()))
        lep.cut=tmpExp
        if(tmpExp.find("\'")>-1):
            lep.cut=tmpExp.strip("'")

