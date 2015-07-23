import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('mcversion', 'Spring14dr', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC campaign or data")
options.register('generator', 'pythia6', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC generator")
options.register('mcWeight', 1.0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "MC sample event weight")
options.register('lepton', 'muon', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Lepton+jets channel")
options.register('metcl', 1, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "MET correction level")

options.register('scaleType', 'abs', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "JES scale type")
options.register('jessource', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Uncertainty source for JES variation for source:up/down")
options.register('flavor', 'bottom', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Jet flavor for flavor:up/down")
options.register('lJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "JES")
options.register('bJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "bJES")
options.register('resolutionSigma', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "JER sigma")
options.register('mesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Muon energy scale shift in sigma")
options.register('eesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Electron energy scale shift in sigma")
options.register('uncFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Unclustered energy factor")

options.register('csvm', 0.679, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "CSVM working point")
options.register('nbjets', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Minimum number of bjets")

options.register('brCorrection', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Do BR correction (MadGraph)")
options.register('bSFNewRecipe', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Use new b-tag SF recipe")
	
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
                              skipEvents = cms.untracked.uint32(0),
                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
                              inputCommands=cms.untracked.vstring(
                                      'keep *',
                                      'drop LHERunInfoProduct_*_*_*',
                                      'drop GenLumiInfoProduct_*_*_*'
                              )
)
if (options.mcversion == "Spring14dr" or options.mcversion == "genLevel"):
  readFiles.extend( [
    '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F4BE4558-BE42-E311-8359-7845C4F9321B.root',
     #'/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/74171050-BF42-E311-9AE0-7845C4FC3A1C.root',
     #'/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0E010BFF-1D43-E311-8CBE-848F69FD287A.root',
     #'/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7609E0BD-1E43-E311-8A6C-00A0D1EE9274.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/00120F7A-84F5-E311-9FBE-002618943910.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/00130EFD-E5F2-E311-9C79-002618943919.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003B7873-00F3-E311-8F81-0025905A48F2.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003CEBCE-21F3-E311-917C-0025905A606A.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/004174FD-9EF5-E311-A561-002354EF3BE4.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/0042B5A2-C1F5-E311-A05D-002618FDA208.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/005D7C0A-11F3-E311-A580-002618943869.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/00673FDD-DCF5-E311-B373-0025905964C2.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/0088B11B-B6F5-E311-BB22-002590596484.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/00972CCE-E5F2-E311-82D9-00259059391E.root',
    #'/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/328C1E89-F9D3-E311-89AA-0017A4770C1C.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00222AFA-04F2-E311-8EE6-0026189438D8.root',
    #'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/4CD5BC5E-CEF5-E311-8A76-0025905A60FE.root',
    #'file:/nfs/dust/cms/user/mseidel/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_Spring14dr_PU20bx25_POSTLS170_V5-v1_AODSIM.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1000_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1002_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1004_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1006_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1008_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_100_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1010_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1012_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1014_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1016_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1018_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1020_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1022_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1024_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1026_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1028_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_102_fastreco_FASTSIM_HLT_PU.root',
    
    #'/store/user/mseidel/TT_Lund_8TeV-sherpa2/job_997_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Lund_8TeV-sherpa2/job_999_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Lund_8TeV-sherpa2/job_99_fastreco_FASTSIM_HLT_PU.root',
    #'/store/user/mseidel/TT_Lund_8TeV-sherpa2/job_9_fastreco_FASTSIM_HLT_PU.root',
  ] )

if (options.mcversion == "Summer12" or options.mcversion == "genLevel"):
  if (options.generator == "pythia8"):
    readFiles.extend( ['/store/user/mseidel/TT_8TeV-amcatnlo-pythia8/job_350_TT_FSIM.root'] )
  elif (options.generator == "sherpa"):
    readFiles.extend( ['/store/user/mseidel/TT_Cluster_8TeV-sherpa2/job_1014_fastreco_FASTSIM_HLT_PU.root'] )
  elif (options.generator == "herwigpp"):
    readFiles.extend( ['/store/user/mseidel/TT_8TeV-amcatnlo-herwigpp/job_35_TT_FSIM.root'] )


#elif (options.mcversion == "Summer12"):
  #readFiles.extend( [
         #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/FC5CECAE-8B14-E211-8578-0025B3E0652A.root'
  #] )
#elif (options.mcversion == "Summer12W0Jets"):
  #readFiles.extend( [
         #'/store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v2/0004/FE9FA8F7-2BF3-E111-A34E-001E672CC1E7.root',
         #'/store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v2/0004/FE7A6C35-FAF2-E111-97BC-D8D385FF4ABA.root',
         #'/store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v2/0004/FE777AEF-CBF2-E111-AD95-001E67398223.root',
  #] )
#elif (data and options.lepton == "muon"):
  #readFiles.extend( [
         #'/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/30000/FEDCB8E2-5270-E211-8FD6-00266CFFBC38.root'
  #] )
#elif (data and options.lepton == "electron"):
  #readFiles.extend( [
         #'/store/data/Run2012A/SingleElectron/AOD/22Jan2013-v1/30002/72352080-9372-E211-B431-00266CFFA1AC.root'
  #] )
secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.source.eventsToProcess = cms.untracked.VEventRange('1:72453536')

## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

## configure geometry & conditions
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if os.getenv('CMSSW_VERSION').startswith('CMSSW_7_1_'):
  process.GlobalTag.globaltag = cms.string('POSTLS170_V6::All')
#elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  #process.GlobalTag.globaltag = cms.string('START42_V17::All')
#elif os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
  #process.GlobalTag.globaltag = cms.string('START53_V27::All')
#if data:
  #process.GlobalTag.globaltag = cms.string('FT_53_V21_AN5::All')

## generator filters
#if (options.mcversion == 'Spring14dr'):
if (options.generator == 'w0jets'):
    process.load("TopAnalysis.TopFilter.filters.GeneratorWNJetsFilter_cfi")
    process.filterWNJets.NJet = 0

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")

## redefine veto jets to be sure it is also replaced when running on PF
#from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJetsPF
process.vetoJets.src = "goodJetsPF30"
process.vetoJets.cut = ''

#process.load("RecoJets.Configuration.GenJetParticles_cff")
#process.load("RecoJets.JetProducers.ak4GenJets_cfi")
#process.load("RecoJets.JetAssociationProducers.kt4JTA_cff")


## std sequence for PAT
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

if not data:
    ## configure JetEnergyScale tool
    process.load("TopAnalysis.TopUtils.JetEnergyScale_cff")
    from TopAnalysis.TopUtils.JetEnergyScale_cff import *

    scaledJetEnergy.scaleType    = options.scaleType
    scaledJetEnergy.JECUncSrcFile= "TopAnalysis/TopUtils/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt" #Uncertainty sources
    scaledJetEnergy.sourceName   = options.jessource
    scaledJetEnergy.flavor       = options.flavor
    scaledJetEnergy.scaleFactor  = options.lJesFactor
    scaledJetEnergy.scaleFactorB = options.bJesFactor

    if options.mcversion == "genLevel":
        resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    else:
        resolutionNominal = [1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056]
    resolutionUnc     = [0.026, 0.028, 0.029, 0.046, 0.062, 0.063, 0.191]

    scaledJetEnergy.resolutionFactors   = []
    for nom,unc in zip(resolutionNominal,resolutionUnc):
      scaledJetEnergy.resolutionFactors.append(nom+options.resolutionSigma*unc)

    scaledJetEnergy.resolutionEtaRanges   = [0.0,0.5, 0.5,1.1, 1.1,1.7, 1.7,2.3, 2.3,2.8, 2.8,3.2, 3.2,-1.]

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
    process.vertexSelectedMuons2012.src   = "scaledMuonEnergy:selectedPatMuons"

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
if os.getenv('CMSSW_VERSION').startswith('CMSSW_7_1_'):
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", throw=True)
#else:
    #process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu24_eta2p1_v*"], throw=True)
#if(options.mcversion=="Summer11"):
    #process.hltFilter.HLTPaths=["HLT_IsoMu24_v*"]
#if(options.mcversion=="Summer12"):
    #process.hltFilter.HLTPaths=["HLT_IsoMu24_eta2p1_v*"]
#if data:
    #process.hltFilter.HLTPaths=["HLT_IsoMu24_eta2p1_v*"]

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
    if (options.mcversion=='genLevel'):
        setForAllTtSemiLepHypotheses(process, "leps", 'goodElectronsEJ')
    else:
        useElectronsForAllTtSemiLepHypotheses(process, 'goodElectronsEJ')
setForAllTtSemiLepHypotheses(process, "jets", "goodJetsPF30")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
setForAllTtSemiLepHypotheses(process, "mets", "scaledMET:scaledMETs")
setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
#if data:
    #setForAllTtSemiLepHypotheses(process, "mets", "patMETs")
    #setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
process.ttSemiLepJetPartonMatch.maxNJets  = -1

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "combinedSecondaryVertexBJetTags"
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = options.csvm
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = options.csvm
#process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 1.0
#process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 3.0
process.hitFitTtSemiLepEventHypothesis.useBTagging       = True

addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )
#if data: removeTtSemiLepHypGenMatch(process)
if (data): removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import analyzeHypothesis
process.analyzeHitFit = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key")
from TopMass.TopEventTree.JetEventAnalyzer_cfi import analyzeJets
process.analyzeJets = analyzeJets.clone(jets = "goodJetsPF20")
from TopMass.TopEventTree.WeightEventAnalyzer_cfi import analyzeWeights
process.analyzeWeights = analyzeWeights.clone(

                                              mcWeight        = options.mcWeight,
                                              puWeightSrc     = cms.InputTag("eventWeightPUsysNo"  , "eventWeightPU"),
                                              puWeightUpSrc   = cms.InputTag("eventWeightPUsysUp"  , "eventWeightPUUp"),
                                              puWeightDownSrc = cms.InputTag("eventWeightPUsysDown", "eventWeightPUDown"),
                                              savePDFWeights = True,
                                              brCorrection   = options.brCorrection
                                             )

process.ttSemiLepHypGenMatch.useBReg = cms.bool(False)
process.ttSemiLepHypMVADisc.useBReg = cms.bool(False)


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

    if (options.generator == "rd"):
        process.eventWeightPU.MCSampleHistoName        = ""
        process.eventWeightPU.DataHistoName            = ""

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
    if options.bSFNewRecipe:
        process.bTagSFEventWeight.newRecipe= True
        process.bTagSFEventWeight.maxJets  = 4
    process.bTagSFEventWeight.sysVar   = "" # bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown possible;
    process.bTagSFEventWeight.filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiency2012.root" 
    process.bTagSFEventWeight.verbose  = 0

    process.bTagSFEventWeightBTagSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
    process.bTagSFEventWeightBTagSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
    process.bTagSFEventWeightMisTagSFUp   = process.bTagSFEventWeight.clone(sysVar = "misTagSFUp")
    process.bTagSFEventWeightMisTagSFDown = process.bTagSFEventWeight.clone(sysVar = "misTagSFDown")


    ## ---
    ##    MC B-JES reweighting
    ## ---
    process.load("TopAnalysis.TopUtils.EventWeightBJES_cfi")

    process.bJESEventWeightFNuUp    = process.EventWeightBJES.clone(
        nuDecayFractionTarget = 0.268
    )
    process.bJESEventWeightFNuDown  = process.EventWeightBJES.clone(
        nuDecayFractionTarget = 0.239
    )
    process.bJESEventWeightFrag     = process.EventWeightBJES.clone(
        fragTargetFile = "TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star_rbLEP.root"
    )
    process.bJESEventWeightFragHard = process.EventWeightBJES.clone(
        fragTargetFile = "TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star_rbLEPhard.root"
    )
    process.bJESEventWeightFragSoft = process.EventWeightBJES.clone(
        fragTargetFile = "TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star_rbLEPsoft.root"
    )

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

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
   maxEventsToPrint = cms.untracked.int32(20),
   printVertex = cms.untracked.bool(False),
   src = cms.InputTag("genParticles")
)


## end path   
#if data:
    #process.path = cms.Path(
                            #process.hltFilter *
                            #process.semiLeptonicSelection *
                            #process.tightBottomSSVPFJets *
                            #process.tightBottomCSVPFJets *
                            #process.semiLeptonicEvents *
                            #process.makeTtSemiLepEvent *
                            #process.analyzeHitFit *
                            #process.analyzeJets *
                            #process.analyzeWeights
                            #)
#else:
if (options.mcversion == "genLevel"):
    delattr(process,"selectedPatJets")
    delattr(process,"goodJetsPF30")
    delattr(process,"goodJetsPF20")
    process.load("TopAnalysis.TopUtils.convertGenToPatJets_cff")
    process.goodJetsPF30 = process.mySelectedPatJets.clone(src = cms.InputTag('scaledJetEnergy','selectedPatJets'),
                                                           cut = ''
                                                          )
    process.goodJetsPF20 = process.mySelectedPatJets.clone(src = cms.InputTag('scaledJetEnergy','selectedPatJets'),
                                                           cut = ''
                                                          )
    #delattr(process,"scaledJetEnergy")
    #process.scaledJetEnergy = process.mySelectedPatJets.clone(src = 'selectedPatJets',
    #                                                       cut = ''
    #                                                      )

    delattr(process,"patElectrons")
    delattr(process,"patMuons")
    delattr(process,"looseElectrons")
    delattr(process,"looseMuons")
    delattr(process,"tightElectronsEJ")
    delattr(process,"tightMuons")
    delattr(process,"patMETs")
    delattr(process,"scaledMET")
    setForAllTtSemiLepHypotheses(process, "mets", "scaledMET")
    process.load("TopAnalysis.TopUtils.convertGenToPatLepton_cff")

    delattr(process,"goodElectronsEJ")
    process.goodElectronsEJ = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("tightElectronsEJ"),
        cut = cms.string('')
    )

    process.convElecRejection.src = 'tightElectronsEJ'

    process.dump = cms.EDAnalyzer('EventContentAnalyzer')

    process.path = cms.Path(#process.printTree*
                            process.convertGenToPatJets *
                            process.scaledJetEnergy *
                            process.convertGenToPatLeptons*
                            #process.dump*
                            process.goodJetsPF20*
                            process.goodJetsPF30*
                            process.tightLeadingPFJets *
                            process.tightBottomCSVPFJets *
                            process.semiLeptonicEvents *
                            process.makeGenEvt *
                            process.analyzeHitFit *
                            process.analyzeJets *
                            process.analyzeWeights
                            )

    #TEST
    #process.ttSemiLepJetPartonMatch.maxDist = 0.5
    #process.ttSemiLepJetPartonMatch.algorithm = "totalMinDist"

    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    #massSearchReplaceAnyInputTag(process.path, cms.InputTag('scaledJetEnergy','selectedPatJets'), 'scaledJetEnergy')
    massSearchReplaceAnyInputTag(process.path, cms.InputTag('scaledMET','scaledMETs'), 'scaledMET')

else:
    process.path = cms.Path(#process.printTree*
                            #process.hltFilter *
                            process.scaledJetEnergy *
                            process.scaledElectronEnergy *
                            process.scaledMuonEnergy *
                            process.scaledMET *
                            process.semiLeptonicSelection *
                            process.tightBottomSSVPFJets *
                            process.tightBottomCSVPFJets *
                            process.semiLeptonicEvents *
                            process.eventWeightMC *
                            #process.makeEventWeightsPU *
                            #process.bTagSFEventWeight *
                            #process.bTagSFEventWeightBTagSFUp *
                            #process.bTagSFEventWeightBTagSFDown *
                            #process.bTagSFEventWeightMisTagSFUp *
                            #process.bTagSFEventWeightMisTagSFDown *
                            process.bJESEventWeightFNuUp *
                            process.bJESEventWeightFNuDown *
                            process.bJESEventWeightFrag *
                            process.bJESEventWeightFragHard *
                            process.bJESEventWeightFragSoft *
                            #process.effSFMuonEventWeight *
                            process.makeGenEvt *
                            #process.makeTtSemiLepEvent *
                            process.analyzeHitFit *
                            process.analyzeJets *
                            process.analyzeWeights
                            )


#delattr(process,"tightBottomJets")
#delattr(process,"tightLeadingJets")

#if (options.mcversion == "Spring14dr"): process.path.remove(process.eventWeightMC)

process.path.remove(process.centralJets)
process.path.remove(process.reliableJets)
process.path.remove(process.goodJets)
process.path.remove(process.trackCountingHighPurBJets)
process.path.remove(process.trackCountingHighEffBJets)
process.path.remove(process.tightLeadingJets)
process.path.remove(process.tightBottomJets)
#process.path.remove(process.)
#process.path.remove(process.unconvTightElectronsEJ)
#process.path.remove(process.goodElectronsEJ)
#process.path.remove(process.looseElectronsEJ)
#process.path.remove(process.tightElectronsEJ)

#if (options.mcversion == "Sherpa12"): process.path.remove(process.eventWeightMC)


## switch to from muon to electron collections
if (options.lepton=="electron"):
    process.analyzeHitFit.lepton = 11
    process.TFileService.fileName = "analyzeTop.root"
    # adpat trigger
    if (options.mcversion == "Fall11"):
        process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
    elif (options.mcversion == "Summer12"):
        process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
    elif (options.mcversion == "Spring14dr"):
        process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
    elif data:
        process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]

    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    pathnames = process.paths_().keys()
    for pathname in pathnames:
      # replace muon selection
      getattr(process, pathname).replace(process.muonSelection, process.electronSelection)
      getattr(process, pathname).remove(process.secondMuonVeto)
      getattr(process, pathname).remove(process.electronVeto)
      ## replace effSF
      if not data: getattr(process, pathname).replace(process.effSFMuonEventWeight, process.effSFElectronEventWeight)
      # replace muon by electron in (remaining) kinfit analyzers
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'tightMuons', 'goodElectronsEJ')
      massSearchReplaceAnyInputTag(getattr(process, pathname), 'effSFMuonEventWeight', 'effSFElectronEventWeight')
      if (options.generator == "sherpa" or options.generator == "herwigpp"):
         #getattr(process, pathname).remove(process.genEvt)
         #getattr(process, pathname).remove(process.decaySubset)
         #getattr(process, pathname).remove(process.initSubset)
         #getattr(process, pathname).remove(process.makeGenEvt)
         getattr(process, pathname).remove(process.eventWeightMC)

if not options.mcversion == "genLevel":
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
        'electronIDs': ['CiC','MVA'],
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
        'addMETSignificance': False,
        }
    if data:
        PFoptions['JECEra' ] = 'Winter14_V2_DATA'
        PFoptions['JECFile'] = '../data/Winter14_DATA_V2PT.db'
        if os.getenv('GC_CONF'):
            print "Running with GC, resetting address of JECFile!"
            PFoptions['JECFile'] = '../src/TopMass/Configuration/data/Winter14_DATA_V2PT.db'

    #if options.mcversion == "Summer12RD":
    if (options.generator == "rd"):
        PFoptions['JECEra' ] = 'Winter14_V1_MC'
        PFoptions['JECFile'] = '../data/Winter14_V1_MC.db'
        if os.getenv('GC_CONF'):
            print "Running with GC, resetting address of JECFile!"
            PFoptions['JECFile'] = '../src/TopMass/Configuration/data/Winter14_V1_MC.db'

    prependPF2PATSequence(process, options = PFoptions)

    ## b-tagging for new jets
    from RecoBTag.Configuration.RecoBTag_cff import impactParameterTagInfos,secondaryVertexTagInfos,combinedSecondaryVertex,combinedSecondaryVertexBJetTags
    #ak4PFCHS
    process.ak4PFCHSImpactParameterTagInfos = impactParameterTagInfos.clone()
    process.ak4PFCHSImpactParameterTagInfos.jetTracks = "jetTracksAssociatorAtVertex"
    process.ak4PFCHSSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
    process.ak4PFCHSSecondaryVertexTagInfos.trackIPTagInfos = "ak4PFCHSImpactParameterTagInfos"
    process.ak4PFCHSStandardCombinedSecondaryVertex = combinedSecondaryVertex.clone()
    process.combinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
    process.combinedSecondaryVertexBJetTags.jetTagComputer = cms.string('ak4PFCHSStandardCombinedSecondaryVertex')
    process.combinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("ak4PFCHSImpactParameterTagInfos"), cms.InputTag("ak4PFCHSSecondaryVertexTagInfos") )
    process.ak4PFCHSJetBtaggingSV = cms.Sequence(process.ak4PFCHSImpactParameterTagInfos *
                                                process.ak4PFCHSSecondaryVertexTagInfos *
                                                process.combinedSecondaryVertexBJetTags
                                                )
    process.ak4PFCHSJetsBtag = cms.Sequence(process.jetTracksAssociatorAtVertex *
                                            process.ak4PFCHSJetBtaggingSV
                                        )
    
    process.patJets.discriminatorSources.remove(cms.InputTag("jetBProbabilityBJetTags"))
    process.patJets.discriminatorSources.remove(cms.InputTag("jetProbabilityBJetTags"))
    process.patJets.discriminatorSources.remove(cms.InputTag("trackCountingHighPurBJetTags"))
    process.patJets.discriminatorSources.remove(cms.InputTag("trackCountingHighEffBJetTags"))
    process.patJets.discriminatorSources.remove(cms.InputTag("simpleSecondaryVertexHighEffBJetTags"))
    process.patJets.discriminatorSources.remove(cms.InputTag("simpleSecondaryVertexHighPurBJetTags"))
    
    process.patJets.addBTagInfo=True
    
    
    
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
        #if (options.mcversion == 'Summer12W0Jets'):
        if (options.generator == 'w0jets'):
            getattr(process, pathname).insert(0,process.filterWNJets)
        
        getattr(process, pathname).insert(0,process.ak4PFCHSJetsBtag)
        getattr(process, pathname).remove(process.goodOfflinePrimaryVertices)
        getattr(process, pathname).insert(0,process.goodOfflinePrimaryVertices)
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
    
