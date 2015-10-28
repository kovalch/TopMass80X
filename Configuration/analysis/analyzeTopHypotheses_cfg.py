import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('mcversion', 'RunIISpring15DR', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC campaign or data")
options.register('runOnMiniAOD'    , True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if run on miniAOD or AOD input' )
options.register('useElecEAIsoCorr', True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if EA (rho) or Delta beta corrections are used for electron isolation is used' )
options.register('useCalibElec'    , False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if electron re-calibration using regression energies is used' )

options.register('generator', 'pythia8', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC generator")
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

options.register('csvm', 0.814, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "CSVM working point")
options.register('nbjets', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Minimum number of bjets")

options.register('brCorrection', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Do BR correction (MadGraph)")
options.register('bSFNewRecipe', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Use new b-tag SF recipe")

options.register('addTriggerMatch' , True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if trigger objects are matched to signal muons' )

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

if (options.mcversion == "data"):
	data = True
	options.register('runOnMC', False , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if run on MC or real data' )
else:
	data = False
	options.register('runOnMC', True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if run on MC or real data' )
## top mass measurement
process = cms.Process("topMass")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.categories.append('TtSemiLeptonicEvent')
process.MessageLogger.cerr.TtSemiLeptonicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *

if data:
    inputFiles = ['/store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/50000/96696580-736F-E511-922D-0025905A60FE.root']
else:
    inputFiles = ['/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/0071E654-216E-E511-9CD4-002590596484.root']


### Selection steps
# If a step is switched off here, its results will still be available in the corresponding TriggerResults of this process.

# Event filter
# This parameter defines the level, at which events will be filtered for output.
# The available levels (paths) are (ordered!):
# 0a. pTrigger
# 0b. pEventCleaning
# 0c. pGoodVertex
# 1.  pSignalMuon
# 2.  pLooseMuonVeto
# 3.  pElectronVeto
# 4a. p1Jet
# 4b. p2Jets
# 4c. p3Jets
# 5.  p4Jets
# 6.  pBTags
# Each level includes the former ones, but also the corresponding stand-alone selection steps are available, adding a
# 'StandAlone' after the prefix 'p' (e.g. 'pLooseMuonVeto' --> 'pStandAloneLooseMuonVeto').
# All corresponding flags are available in the TriggerResults collection produced by this process later.
selectEvents = 'pGoodVertex'

# Step 0
#triggerSelectionData = ''
#triggerSelectionMC   = ''

# Step 1
#muonCut       = ''
#signalMuonCut = ''
#muonVertexMaxDZ = 0.5

# Step 2

# Step 3
useElecEAIsoCorr = options.useElecEAIsoCorr
useCalibElec     = options.useCalibElec
#electronGsfCut   = ''
#electronCalibCut = ''
electronCut = electronGsfCut

# Step 4

#jetCut = ''
# Step4a
#veryTightJetCut = ''
# Step4b
#tightJetCut     = ''
# Step4c
#looseJetCut     = ''

# Step 5
#veryLooseJetCut = ''

# Step 6
bTagSrc = 'selectedJets'
bTagCut = 'bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.814'
minBTags = 2

# TriggerMatching
addTriggerMatch = options.addTriggerMatch
#triggerObjectSelectionData = 'type("TriggerMuon") && ( path("%s") )'%( triggerSelectionData )
#triggerObjectSelectionMC   = 'type("TriggerMuon") && ( path("%s") )'%( triggerSelectionMC )


### Input

runOnMC      = options.runOnMC
runOnMiniAOD = options.runOnMiniAOD

# maximum number of events
maxEvents = options.maxEvents

### Conditions

# GlobalTags
#globalTagMC   = '74X_mcRun2_asymptotic_v2'
globalTagData = '74X_dataRun2_v4'
globalTagMC   = 'DEFAULT'
#globalTagData = 'DEFAULT'
usePrivateSQlite=False #use external JECs (sqlite file)

### Output

# output file
outputFile = 'patRefSel_muJets.root'

# event frequency of Fwk report
fwkReportEvery = max( 1000, int( maxEvents / 100 ) )

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True

### ======================================================================== ###
###                                                                          ###
###                              End of constants                            ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


triggerSelection       = triggerSelectionData
triggerObjectSelection = triggerObjectSelectionData
if runOnMC:
  triggerSelection       = triggerSelectionMC
  triggerObjectSelection = triggerObjectSelectionMC


###
### Basic configuration
###
from FWCore.MessageService.MessageLogger_cfi import *
process.options = cms.untracked.PSet(
  wantSummary      = cms.untracked.bool( wantSummary ),
  allowUnscheduled = cms.untracked.bool( True )
)

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery

if runOnMC:
  if globalTagMC == 'DEFAULT':
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
    process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:run2_mc' )
  else:
    process.GlobalTag.globaltag = globalTagMC
else:
  if globalTagData == 'DEFAULT':
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
    process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:run2_data' )
  else:
    process.GlobalTag.globaltag = globalTagData

#override JEC
if usePrivateSQlite:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    if runOnMC:
        era="Summer15_50nsV4_MC"
    else:
        era="Summer15_50nsV4_DATA"
    dBFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/PatAlgos/test/"+era+".db")
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string( "sqlite_file://"+dBFile ),
                               toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
                ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
                ),
            )
                               )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

#redo jets and MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
if runOnMiniAOD:
    runMetCorAndUncFromMiniAOD( process, isData=(runOnMC != True))
else:
    runMETCorrectionsAndUncertainties( process, isData=(runOnMC != True) )
#add mnodules to redo the b-tagging on-the-fly using unscheduled mode
from PhysicsTools.PatAlgos.tools.jetTools import *
## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

jetCorrectionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
if runOnMC == False:
    jetCorrectionLevels.append('L2L3Residual')

from PhysicsTools.PatAlgos.tools.jetTools import *
## call addJetCollection to get get b tagging modules prepared
addJetCollection(
    process,
    labelName = '',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', jetCorrectionLevels, 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #getJetMCFlavour = runOnMC, //TODO does not work at the moment
    getJetMCFlavour = False,
    algo = 'AK',
    rParam = 0.4
)

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
if runOnMiniAOD:
    adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))


#process.patJets.addBTagInfo = cms.bool(True)
#if runOnMC == False:
#   process.patJetCorrFactors.levels.append('L2L3Residual')
###
### Input configuration
###

if len( inputFiles ) == 0:
  if runOnMiniAOD:
    if runOnMC:
      from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValTTbarPileUpMINIAODSIM
      inputFiles = filesRelValTTbarPileUpMINIAODSIM
    else:
      from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValSingleMuMINIAOD
      inputFiles = filesRelValSingleMuMINIAOD
  else:
    if runOnMC:
      from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
      inputFiles = filesRelValProdTTbarAODSIM
    else:
      from PhysicsTools.PatAlgos.patInputFiles_cff import filesSingleMuRECO # not available at CERN
      inputFiles = filesSingleMuRECO

process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.source.fileNames = inputFiles
process.maxEvents.input  = maxEvents


###
### PAT configuration
###

if not runOnMiniAOD:
  process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )



###
### Output configuration
###
#
# process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
# # output file name
# process.out.fileName = outputFile
# from TopQuarkAnalysis.Configuration.patRefSel_eventContent_cff import refMuJets_eventContent
# process.out.outputCommands += refMuJets_eventContent
# if runOnMiniAOD:
#   from TopQuarkAnalysis.Configuration.patRefSel_eventContent_cff import miniAod_eventContent
#   process.out.outputCommands += miniAod_eventContent
# else:
#   from TopQuarkAnalysis.Configuration.patRefSel_eventContent_cff import aod_eventContent
#   process.out.outputCommands += aod_eventContent
# # clear event selection
# process.out.SelectEvents.SelectEvents = cms.vstring( selectEvents )


## generator filters
#if (options.mcversion == 'Spring14dr'):
if (options.generator == 'w0jets'):
    process.load("TopAnalysis.TopFilter.filters.GeneratorWNJetsFilter_cfi")
    process.filterWNJets.NJet = 0

###
### Selection configuration
###

# Individual steps

# Step 0

from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.triggerSelection = triggerResults.clone( triggerConditions = [ triggerSelection ] )
process.sStandAloneTrigger = cms.Sequence( process.triggerSelection
                                         )
process.pStandAloneTrigger = cms.Path( process.sStandAloneTrigger )

process.load('TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff')
process.sStandAloneEventCleaning = cms.Sequence()
if runOnMiniAOD:
  process.sStandAloneEventCleaning += process.eventCleaningMiniAOD
  if runOnMC:
    process.sStandAloneEventCleaning += process.eventCleaningMiniAODMC
  else:
    process.sStandAloneEventCleaning += process.eventCleaningMiniAODData
    process.metFiltersMiniAOD.TriggerResultsTag = cms.InputTag("TriggerResults","","RECO")
    #HS disable event cleaning at the moment
    process.sStandAloneEventCleaning = cms.Sequence()
else:
  process.sStandAloneEventCleaning += process.eventCleaning
  if runOnMC:
    process.sStandAloneEventCleaning += process.eventCleaningMC
  else:
    process.sStandAloneEventCleaning += process.eventCleaningData
process.pStandAloneEventCleaning = cms.Path( process.sStandAloneEventCleaning )

from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
process.goodOfflinePrimaryVertices = goodOfflinePrimaryVertices.clone( filter = True )
if runOnMiniAOD:
  process.goodOfflinePrimaryVertices.src = 'offlineSlimmedPrimaryVertices'
process.sStandAloneGoodVertex = cms.Sequence( process.goodOfflinePrimaryVertices
                                            )
process.pStandAloneGoodVertex = cms.Path( process.sStandAloneGoodVertex )


# Step 1

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import selectedMuons, preSignalMuons, signalMuons, standAloneSignalMuonFilter
process.selectedMuons = selectedMuons.clone( cut = muonCut )
if runOnMiniAOD:
  process.selectedMuons.src = 'slimmedMuons'
process.preSignalMuons = preSignalMuons.clone( cut = signalMuonCut )
process.signalMuons = signalMuons.clone( maxDZ = muonVertexMaxDZ )
if runOnMiniAOD:
  process.signalMuons.vertexSource = 'offlineSlimmedPrimaryVertices'
process.standAloneSignalMuonFilter = standAloneSignalMuonFilter.clone()
process.sStandAloneSignalMuon = cms.Sequence( process.standAloneSignalMuonFilter )
process.pStandAloneSignalMuon = cms.Path( process.sStandAloneSignalMuon )

# Step 2

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import standAloneLooseMuonVetoFilter
process.standAloneLooseMuonVetoFilter = standAloneLooseMuonVetoFilter.clone()
process.sStandAloneLooseMuonVeto = cms.Sequence( process.standAloneLooseMuonVetoFilter )
process.pStandAloneLooseMuonVeto = cms.Path( process.sStandAloneLooseMuonVeto )

# Step 3

if not runOnMiniAOD:
  from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection
  electron_ids = [ 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_CSA14_50ns_V1_cff'
                 , 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_CSA14_PU20bx25_V0_cff'
                  ]
  switchOnVIDElectronIdProducer( process )
  process.electronIDValueMapProducer.ebReducedRecHitCollection = cms.InputTag( 'reducedEcalRecHitsEB' )
  process.electronIDValueMapProducer.eeReducedRecHitCollection = cms.InputTag( 'reducedEcalRecHitsEE' )
  process.electronIDValueMapProducer.esReducedRecHitCollection = cms.InputTag( 'reducedEcalRecHitsES' )
  for idmod in electron_ids:
    setupAllVIDIdsInModule( process, idmod, setupVIDElectronSelection )

if useElecEAIsoCorr:
  from EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi import elPFIsoValueEA03
  if runOnMiniAOD:
    process.patElPFIsoValueEA03 = elPFIsoValueEA03.clone( gsfElectrons = ''
                                                        , pfElectrons  = ''
                                                        , patElectrons = cms.InputTag( 'slimmedElectrons' )
                                                        , rhoIso       = cms.InputTag( 'fixedGridRhoFastjetAll' )
                                                        )
    from EgammaAnalysis.ElectronTools.patElectronEAIsoCorrectionProducer_cfi import patElectronEAIso03CorrectionProducer
    process.electronsWithEA03Iso = patElectronEAIso03CorrectionProducer.clone( patElectrons  = 'slimmedElectrons'
                                                                             , eaIsolator    = 'patElPFIsoValueEA03'
                                                                             )
  else:
    process.elPFIsoValueEA03 = elPFIsoValueEA03.clone( gsfElectrons = 'gedGsfElectrons'
                                                     , pfElectrons  = ''
                                                     , rhoIso       = cms.InputTag( 'fixedGridRhoFastjetAll' )
                                                     )
    process.patElectrons.isolationValues.user = cms.VInputTag( cms.InputTag( 'elPFIsoValueEA03' ) )
else:
  electronGsfCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )
  electronCalibCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )
  electronCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )

if useCalibElec:
  from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import electronsWithRegression, calibratedElectrons
  process.electronsWithRegression = electronsWithRegression.clone()
  if runOnMiniAOD:
    if useElecEAIsoCorr:
      process.electronsWithRegression.inputElectronsTag = 'electronsWithEA03Iso'
    else:
      process.electronsWithRegression.inputElectronsTag = 'slimmedElectrons'
    process.electronsWithRegression.vertexCollection  = 'offlineSlimmedPrimaryVertices'
  process.calibratedElectrons = calibratedElectrons.clone( isMC = runOnMC )
  if runOnMC:
    process.calibratedElectrons.inputDataset = 'Summer12_LegacyPaper' # FIXME: Update as soon as available
  else:
    process.calibratedElectrons.inputDataset = '22Jan2013ReReco' # FIXME: Update as soon as available
  process.RandomNumberGeneratorService = cms.Service( "RandomNumberGeneratorService"
                                                    , calibratedElectrons = cms.PSet( initialSeed = cms.untracked.uint32( 1 )
                                                                                    , engineName  = cms.untracked.string('TRandom3')
                                                                                    )
                                                    )
  electronCut = electronCalibCut



from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import selectedElectrons, standAloneElectronVetoFilter
process.selectedElectrons = selectedElectrons.clone( cut = electronCut )
if useCalibElec:
  process.selectedElectrons.src = 'calibratedElectrons'
elif useElecEAIsoCorr and runOnMiniAOD:
  process.selectedElectrons.src = 'electronsWithEA03Iso'
elif runOnMiniAOD:
  process.selectedElectrons.src = 'slimmedElectrons'

process.standAloneElectronVetoFilter = standAloneElectronVetoFilter.clone()
process.sStandAloneElectronVeto = cms.Sequence( process.standAloneElectronVetoFilter )
process.pStandAloneElectronVeto = cms.Path( process.sStandAloneElectronVeto )

# Step 4

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import selectedJets
process.selectedJets = selectedJets.clone( cut = jetCut )
#if runOnMiniAOD:
#  process.selectedJets.src = 'slimmedJets'
process.selectedJets.src = 'patJets'

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import signalVeryTightJets, standAloneSignalVeryTightJetsFilter
process.signalVeryTightJets = signalVeryTightJets.clone( cut = veryTightJetCut )
process.standAloneSignalVeryTightJetsFilter = standAloneSignalVeryTightJetsFilter.clone()
process.sStandAlone1Jet = cms.Sequence( process.standAloneSignalVeryTightJetsFilter )
process.pStandAlone1Jet = cms.Path( process.sStandAlone1Jet )

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import signalTightJets, standAloneSignalTightJetsFilter
process.signalTightJets = signalTightJets.clone( cut = tightJetCut )
process.standAloneSignalTightJetsFilter = standAloneSignalTightJetsFilter.clone()
process.sStandAlone2Jets = cms.Sequence( process.standAloneSignalTightJetsFilter )
process.pStandAlone2Jets = cms.Path( process.sStandAlone2Jets )

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import signalLooseJets, standAloneSignalLooseJetsFilter
process.signalLooseJets = signalLooseJets.clone( cut = looseJetCut )
process.standAloneSignalLooseJetsFilter = standAloneSignalLooseJetsFilter.clone()
process.sStandAlone3Jets = cms.Sequence( process.standAloneSignalLooseJetsFilter )
process.pStandAlone3Jets = cms.Path( process.sStandAlone3Jets )

# Step 5

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import signalVeryLooseJets, standAloneSignalVeryLooseJetsFilter
process.signalVeryLooseJets = signalVeryLooseJets.clone( cut = veryLooseJetCut )
process.standAloneSignalVeryLooseJetsFilter = standAloneSignalVeryLooseJetsFilter.clone()
process.sStandAlone4Jets = cms.Sequence( process.standAloneSignalVeryLooseJetsFilter )
process.pStandAlone4Jets = cms.Path( process.sStandAlone4Jets )

# Step 6

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import selectedBTagJets, standAloneSignalBTagsFilter
process.selectedBTagJets = selectedBTagJets.clone( src = bTagSrc
                                                 , cut = bTagCut
                                                 )
process.standAloneSignalBTagsFilter = standAloneSignalBTagsFilter.clone( minNumber = minBTags )
process.sStandAloneBTags = cms.Sequence( process.standAloneSignalBTagsFilter )
process.pStandAloneBTags = cms.Path( process.sStandAloneBTags )

# Consecutive steps

process.sTrigger       = cms.Sequence( process.sStandAloneTrigger
                                     )
process.sEventCleaning = cms.Sequence( process.sTrigger
                                     + process.sStandAloneEventCleaning
                                     )
process.sGoodVertex    = cms.Sequence( process.sEventCleaning
                                     + process.sStandAloneGoodVertex
                                     )
process.sSignalMuon    = cms.Sequence( process.sGoodVertex
                                     + process.sStandAloneSignalMuon
                                     )
process.sLooseMuonVeto = cms.Sequence( process.sSignalMuon
                                     + process.sStandAloneLooseMuonVeto
                                     )
process.sElectronVeto  = cms.Sequence( process.sLooseMuonVeto
                                     + process.sStandAloneElectronVeto
                                     )
process.s1Jet          = cms.Sequence( process.sElectronVeto
                                     + process.sStandAlone1Jet
                                     )
process.s2Jets         = cms.Sequence( process.s1Jet
                                     + process.sStandAlone2Jets
                                     )
process.s3Jets         = cms.Sequence( process.s2Jets
                                     + process.sStandAlone2Jets
                                     )
process.s4Jets         = cms.Sequence( process.s3Jets
                                     + process.sStandAlone4Jets
                                     )
process.sBTags         = cms.Sequence( process.s4Jets
                                     + process.sStandAloneBTags
                                     )

process.pTrigger       = cms.Path( process.sTrigger )
process.pEventCleaning = cms.Path( process.sEventCleaning )
process.pGoodVertex    = cms.Path( process.sGoodVertex )
process.pSignalMuon    = cms.Path( process.sSignalMuon )
process.pLooseMuonVeto = cms.Path( process.sLooseMuonVeto )
process.pElectronVeto  = cms.Path( process.sElectronVeto )
process.p1Jet          = cms.Path( process.s1Jet )
process.p2Jets         = cms.Path( process.s2Jets )
process.p3Jets         = cms.Path( process.s2Jets )
process.p4Jets         = cms.Path( process.s4Jets )


process.pBTags         = cms.Path( process.sBTags )

# Trigger matching

if addTriggerMatch:
  from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cff import muonTriggerMatch
  process.muonTriggerMatch = muonTriggerMatch.clone( matchedCuts = triggerObjectSelection )
  if not runOnMiniAOD:
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTriggerMatchEmbedding
    switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'muonTriggerMatch' ] )
  else:
    from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cff import unpackedPatTrigger
    process.selectedTriggerUnpacked = unpackedPatTrigger.clone()
    process.muonTriggerMatch.matched = 'selectedTriggerUnpacked'
    from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cff import signalMuonsTriggerMatch
    process.signalMuonsTriggerMatch = signalMuonsTriggerMatch.clone()
    # process.out.outputCommands += [ 'drop *_signalMuons_*_*'
    #                               , 'keep *_signalMuonsTriggerMatch_*_*'
    #                               ]





# if not data:
#     ## configure JetEnergyScale tool
#     process.load("TopAnalysis.TopUtils.JetEnergyScale_cff")
#     from TopAnalysis.TopUtils.JetEnergyScale_cff import *
#
#     scaledJetEnergy.scaleType    = options.scaleType
#     scaledJetEnergy.JECUncSrcFile= "TopAnalysis/TopUtils/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt" #Uncertainty sources
#     scaledJetEnergy.sourceName   = options.jessource
#     scaledJetEnergy.flavor       = options.flavor
#     scaledJetEnergy.scaleFactor  = options.lJesFactor
#     scaledJetEnergy.scaleFactorB = options.bJesFactor
#
#     if options.mcversion == "genLevel":
#         resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#     else:
#         resolutionNominal = [1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056]
#     resolutionUnc     = [0.026, 0.028, 0.029, 0.046, 0.062, 0.063, 0.191]
#
#     scaledJetEnergy.resolutionFactors   = []
#     for nom,unc in zip(resolutionNominal,resolutionUnc):
#       scaledJetEnergy.resolutionFactors.append(nom+options.resolutionSigma*unc)
#
#     scaledJetEnergy.resolutionEtaRanges   = [0.0,0.5, 0.5,1.1, 1.1,1.7, 1.7,2.3, 2.3,2.8, 2.8,3.2, 3.2,-1.]
#
#     scaledJetEnergy.inputJets    = "selectedJets"
#     scaledJetEnergy.inputMETs    = "patMETs"
#
#     process.signalJets.src = "scaledJetEnergy:selectedJets"
#
#     ## electron shift
#     process.load("TopAnalysis.TopUtils.ElectronEnergyScale_cfi")
#     from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *
#
#     process.scaledElectronEnergy.src      = "selectedElectrons"
#     process.scaledElectronEnergy.mets     = "scaledJetEnergy:patMETs"
#     process.scaledElectronEnergy.shiftBy  = options.eesShift
#     process.signalElectrons.src   = "scaledElectronEnergy:selectedPatElectrons"
#
#     ## muon shift
#     process.load("TopAnalysis.TopUtils.MuonEnergyScale_cfi")
#     from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *
#
#     process.scaledMuonEnergy.src      = "selectedMuons"
#     process.scaledMuonEnergy.mets     = "scaledElectronEnergy:METs"
#     process.scaledMuonEnergy.shiftBy  = options.mesShift
#     process.signalMuons.src   = "scaledMuonEnergy:selectedMuons"
#
#     ## unclustered energy scale
#     process.load("TopAnalysis.TopUtils.UnclusteredMETScale_cfi")
#     from TopAnalysis.TopUtils.UnclusteredMETScale_cfi import *
#
#     process.scaledMET.inputJets       = "scaledJetEnergy:selectedPatJets"
#     process.scaledMET.inputMETs       = "scaledMuonEnergy:METs"
#     process.scaledMET.inputElectrons  = "scaledElectronEnergy:selectedPatElectrons"
#     process.scaledMET.inputMuons      = "scaledMuonEnergy:selectedPatMuons"
#     process.scaledMET.scaleFactor     = options.uncFactor

    ## sequence for ttGenEvent
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")


## sequence for TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

## TRIGGER
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if os.getenv('CMSSW_VERSION').startswith('CMSSW_7_4_'):
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
# process.tightBottomSSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
#                                            cut='bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") > 1.74'
#                                            )
# process.tightBottomCSVPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30',
#                                            cut='bDiscriminator(\"combinedSecondaryVertexBJetTags\") > ' + str(options.csvm)
#                                            )
#
# process.leadingJetSelection.src = 'tightLeadingPFJets'
#
# process.bottomJetSelection.src  = 'tightBottomCSVPFJets'
# process.bottomJetSelection.minNumber = options.nbjets;

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
setForAllTtSemiLepHypotheses(process, "leps", "signalMuons")
if (options.lepton=='electron'):
    if (options.mcversion=='genLevel'):
        setForAllTtSemiLepHypotheses(process, "leps", 'goodElectronsEJ')
    else:
        useElectronsForAllTtSemiLepHypotheses(process, 'signalElectrons')
setForAllTtSemiLepHypotheses(process, "jets", "selectedJets")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
#setForAllTtSemiLepHypotheses(process, "mets", "scaledMET:scaledMETs")
setForAllTtSemiLepHypotheses(process, "mets", "slimmedMETs")

setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
#if data:
    #setForAllTtSemiLepHypotheses(process, "mets", "patMETs")
    #setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")

## change jet-parton matching algorithm
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
process.ttSemiLepJetPartonMatch.maxNJets  = -1

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = "pfCombinedInclusiveSecondaryVertexV2BJetTags"
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
process.analyzeJets = analyzeJets.clone(jets = "selectedJets")
if runOnMiniAOD:
    process.analyzeJets.met = "slimmedMETs"
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


if runOnMiniAOD:
	process.initSubset.src = "prunedGenParticles"
	process.decaySubset.src = "prunedGenParticles"


#process.analyzer         = cms.Sequence(process.makeGenEvt+process.makeTtSemiLepEvent+process.analyzeHitFit+process.analyzeJets+process.analyzeWeights )
process.analyzer         = cms.Sequence(process.analyzeHitFit+process.analyzeJets+process.analyzeWeights )
if data:
    process.pAnal            = cms.Path(process.sElectronVeto+process.analyzer )
else:
    process.pAnal            = cms.Path(process.analyzer )

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
    elif (options.mcversion == "RunIISpring15DR"):
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
    process.bTagSFEventWeight.jets     = "signalTightJets"
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
    if runOnMiniAOD:
        process.EventWeightBJES.genJets      = "slimmedGenJets"
        process.EventWeightBJES.genParticles = "prunedGenParticles"

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
    process.effSFMuonEventWeight.particles=cms.InputTag("signalMuons")
    process.effSFMuonEventWeight.sysVar   = cms.string("")
    process.effSFMuonEventWeight.filename= "TopAnalysis/Configuration/data/MuonEffSF2D2012.root"
    process.effSFMuonEventWeight.verbose=cms.int32(0)
    process.effSFMuonEventWeight.additionalFactor=1.0 ## lepton selection and trigger eff. SF both included in loaded histo
    process.effSFMuonEventWeight.additionalFactorErr=0.01 ## 1.0% sys error to account for selection difference Z - ttbar
    process.effSFMuonEventWeight.shapeDistortionFactor=-1

    process.effSFElectronEventWeight=process.effSFLepton2DEventWeight.clone()
    process.effSFElectronEventWeight.particles=cms.InputTag("signalElectrons")
    process.effSFElectronEventWeight.jets=cms.InputTag("signalTightJets")
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
# if (options.mcversion == "genLevel"):
#     delattr(process,"selectedPatJets")
#     delattr(process,"goodJetsPF30")
#     delattr(process,"goodJetsPF20")
#     process.load("TopAnalysis.TopUtils.convertGenToPatJets_cff")
#     process.goodJetsPF30 = process.mySelectedPatJets.clone(src = cms.InputTag('scaledJetEnergy','selectedPatJets'),
#                                                            cut = ''
#                                                           )
#     process.goodJetsPF20 = process.mySelectedPatJets.clone(src = cms.InputTag('scaledJetEnergy','selectedPatJets'),
#                                                            cut = ''
#                                                           )
#     #delattr(process,"scaledJetEnergy")
#     #process.scaledJetEnergy = process.mySelectedPatJets.clone(src = 'selectedPatJets',
#     #                                                       cut = ''
#     #                                                      )
#
#     delattr(process,"patElectrons")
#     delattr(process,"patMuons")
#     delattr(process,"looseElectrons")
#     delattr(process,"looseMuons")
#     delattr(process,"tightElectronsEJ")
#     delattr(process,"tightMuons")
#     delattr(process,"patMETs")
#     delattr(process,"scaledMET")
#     setForAllTtSemiLepHypotheses(process, "mets", "scaledMET")
#     process.load("TopAnalysis.TopUtils.convertGenToPatLepton_cff")
#
#     delattr(process,"goodElectronsEJ")
#     process.goodElectronsEJ = cms.EDFilter("PATMuonSelector",
#         src = cms.InputTag("tightElectronsEJ"),
#         cut = cms.string('')
#     )
#
#     process.convElecRejection.src = 'tightElectronsEJ'
#
#     process.dump = cms.EDAnalyzer('EventContentAnalyzer')
#
#     process.path = cms.Path(#process.printTree*
#                             process.convertGenToPatJets *
#                             process.scaledJetEnergy *
#                             process.convertGenToPatLeptons*
#                             #process.dump*
#                             process.goodJetsPF20*
#                             process.goodJetsPF30*
#                             process.tightLeadingPFJets *
#                             process.tightBottomCSVPFJets *
#                             process.semiLeptonicEvents *
#                             process.makeGenEvt *
#                             process.analyzeHitFit *
#                             process.analyzeJets *
#                             process.analyzeWeights
#                             )
#
#     #TEST
#     #process.ttSemiLepJetPartonMatch.maxDist = 0.5
#     #process.ttSemiLepJetPartonMatch.algorithm = "totalMinDist"
#
#     from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#     #massSearchReplaceAnyInputTag(process.path, cms.InputTag('scaledJetEnergy','selectedPatJets'), 'scaledJetEnergy')
#     massSearchReplaceAnyInputTag(process.path, cms.InputTag('scaledMET','scaledMETs'), 'scaledMET')
#
# else:
#     process.path = cms.Path(#process.printTree*
#                             #process.hltFilter *
#                             process.scaledJetEnergy *
#                             process.scaledElectronEnergy *
#                             process.scaledMuonEnergy *
#                             process.scaledMET *
#                             process.semiLeptonicSelection *
#                             process.tightBottomSSVPFJets *
#                             process.tightBottomCSVPFJets *
#                             process.semiLeptonicEvents *
#                             process.eventWeightMC *
#                             #process.makeEventWeightsPU *
#                             #process.bTagSFEventWeight *
#                             #process.bTagSFEventWeightBTagSFUp *
#                             #process.bTagSFEventWeightBTagSFDown *
#                             #process.bTagSFEventWeightMisTagSFUp *
#                             #process.bTagSFEventWeightMisTagSFDown *
#                             process.bJESEventWeightFNuUp *
#                             process.bJESEventWeightFNuDown *
#                             process.bJESEventWeightFrag *
#                             process.bJESEventWeightFragHard *
#                             process.bJESEventWeightFragSoft *
#                             #process.effSFMuonEventWeight *
#                             process.makeGenEvt *
#                             #process.makeTtSemiLepEvent *
#                             process.analyzeHitFit *
#                             process.analyzeJets *
#                             process.analyzeWeights
#                             )
#
#
# #delattr(process,"tightBottomJets")
# #delattr(process,"tightLeadingJets")
#
# #if (options.mcversion == "Spring14dr"): process.path.remove(process.eventWeightMC)
#
# process.path.remove(process.centralJets)
# process.path.remove(process.reliableJets)
# process.path.remove(process.goodJets)
# process.path.remove(process.trackCountingHighPurBJets)
# process.path.remove(process.trackCountingHighEffBJets)
# process.path.remove(process.tightLeadingJets)
# process.path.remove(process.tightBottomJets)
# #process.path.remove(process.)
# #process.path.remove(process.unconvTightElectronsEJ)
# #process.path.remove(process.goodElectronsEJ)
# #process.path.remove(process.looseElectronsEJ)
# #process.path.remove(process.tightElectronsEJ)
#
# #if (options.mcversion == "Sherpa12"): process.path.remove(process.eventWeightMC)
#

## switch to from muon to electron collections
# if (options.lepton=="electron"):
#     process.analyzeHitFit.lepton = 11
#     process.TFileService.fileName = "analyzeTop.root"
#     # adpat trigger
#     if (options.mcversion == "Fall11"):
#         process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
#     elif (options.mcversion == "Summer12"):
#         process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
#     elif (options.mcversion == "Spring14dr"):
#         process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
#     elif data:
#         process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
#
#     from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#     pathnames = process.paths_().keys()
#     for pathname in pathnames:
#       # replace muon selection
#       getattr(process, pathname).replace(process.muonSelection, process.electronSelection)
#       getattr(process, pathname).remove(process.secondMuonVeto)
#       getattr(process, pathname).remove(process.electronVeto)
#       ## replace effSF
#       if not data: getattr(process, pathname).replace(process.effSFMuonEventWeight, process.effSFElectronEventWeight)
#       # replace muon by electron in (remaining) kinfit analyzers
#       massSearchReplaceAnyInputTag(getattr(process, pathname), 'tightMuons', 'goodElectronsEJ')
#       massSearchReplaceAnyInputTag(getattr(process, pathname), 'effSFMuonEventWeight', 'effSFElectronEventWeight')
#       if (options.generator == "sherpa" or options.generator == "herwigpp"):
#          #getattr(process, pathname).remove(process.genEvt)
#          #getattr(process, pathname).remove(process.decaySubset)
#          #getattr(process, pathname).remove(process.initSubset)
#          #getattr(process, pathname).remove(process.makeGenEvt)
#          getattr(process, pathname).remove(process.eventWeightMC)

# if not options.mcversion == "genLevel":
#     from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
#     PFoptions = {
#         'runOnMC': not data,
#         'runOnAOD': True,
#         'switchOffEmbedding': False,
#         'addResolutions': True,
#         'resolutionsVersion': 'fall11',
#         'runOnOLDcfg': True,
#         'cutsMuon': 'pt > 10. & abs(eta) < 2.5',
#         'cutsElec': 'et > 20. & abs(eta) < 2.5',
#         'cutsJets': 'pt > 10 & abs(eta) < 5.0',
#         'electronIDs': ['CiC','MVA'],
#         'pfIsoConeMuon': 0.4,
#         'pfIsoConeElec': 0.3,
#         'pfIsoValMuon': 0.2,
#         'pfIsoValElec': 0.15,
#         'doDeltaBetaCorrMuon' : True,
#         'doDeltaBetaCorrElec' : True,
#         'skipIfNoPFMuon': False,
#         'skipIfNoPFElec': False,
#         'addNoCutPFMuon': False,
#         'addNoCutPFElec': False,
#         'noMuonTopProjection': False,
#         'noElecTopProjection': False,
#         'analyzersBeforeMuonIso':cms.Sequence(),
#         'analyzersBeforeElecIso':cms.Sequence(),
#         'excludeElectronsFromWsFromGenJets': True,
#         'METCorrectionLevel': options.metcl,
#         'addMETSignificance': False,
#         }
#     if data:
#         PFoptions['JECEra' ] = 'Winter14_V2_DATA'
#         PFoptions['JECFile'] = '../data/Winter14_DATA_V2PT.db'
#         if os.getenv('GC_CONF'):
#             print "Running with GC, resetting address of JECFile!"
#             PFoptions['JECFile'] = '../src/TopMass/Configuration/data/Winter14_DATA_V2PT.db'
#
#     #if options.mcversion == "Summer12RD":
#     if (options.generator == "rd"):
#         PFoptions['JECEra' ] = 'Winter14_V1_MC'
#         PFoptions['JECFile'] = '../data/Winter14_V1_MC.db'
#         if os.getenv('GC_CONF'):
#             print "Running with GC, resetting address of JECFile!"
#             PFoptions['JECFile'] = '../src/TopMass/Configuration/data/Winter14_V1_MC.db'
#
#     prependPF2PATSequence(process, options = PFoptions)
#
#     ## b-tagging for new jets
#     from RecoBTag.Configuration.RecoBTag_cff import impactParameterTagInfos,secondaryVertexTagInfos,combinedSecondaryVertex,combinedSecondaryVertexBJetTags
#     #ak4PFCHS
#     process.ak4PFCHSImpactParameterTagInfos = impactParameterTagInfos.clone()
#     process.ak4PFCHSImpactParameterTagInfos.jetTracks = "jetTracksAssociatorAtVertex"
#     process.ak4PFCHSSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
#     process.ak4PFCHSSecondaryVertexTagInfos.trackIPTagInfos = "ak4PFCHSImpactParameterTagInfos"
#     process.ak4PFCHSStandardCombinedSecondaryVertex = combinedSecondaryVertex.clone()
#     process.combinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
#     process.combinedSecondaryVertexBJetTags.jetTagComputer = cms.string('ak4PFCHSStandardCombinedSecondaryVertex')
#     process.combinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("ak4PFCHSImpactParameterTagInfos"), cms.InputTag("ak4PFCHSSecondaryVertexTagInfos") )
#     process.ak4PFCHSJetBtaggingSV = cms.Sequence(process.ak4PFCHSImpactParameterTagInfos *
#                                                 process.ak4PFCHSSecondaryVertexTagInfos *
#                                                 process.combinedSecondaryVertexBJetTags
#                                                 )
#     process.ak4PFCHSJetsBtag = cms.Sequence(process.jetTracksAssociatorAtVertex *
#                                             process.ak4PFCHSJetBtaggingSV
#                                         )
#
#     process.patJets.discriminatorSources.remove(cms.InputTag("jetBProbabilityBJetTags"))
#     process.patJets.discriminatorSources.remove(cms.InputTag("jetProbabilityBJetTags"))
#     process.patJets.discriminatorSources.remove(cms.InputTag("trackCountingHighPurBJetTags"))
#     process.patJets.discriminatorSources.remove(cms.InputTag("trackCountingHighEffBJetTags"))
#     process.patJets.discriminatorSources.remove(cms.InputTag("simpleSecondaryVertexHighEffBJetTags"))
#     process.patJets.discriminatorSources.remove(cms.InputTag("simpleSecondaryVertexHighPurBJetTags"))
#
#     process.patJets.addBTagInfo=True
#
#
#
#     ## adaptions (re-aranging of modules) to speed up processing
#     pathnames = process.paths_().keys()
#     for pathname in pathnames:
#         if not data:
#             ## move the ttGenEvent filter to the beginning of the sequence
#             getattr(process, pathname).remove(process.genEvt)
#             getattr(process, pathname).insert(0,process.genEvt)
#             getattr(process, pathname).remove(process.decaySubset)
#             getattr(process, pathname).insert(0,process.decaySubset)
#             getattr(process, pathname).remove(process.initSubset)
#             getattr(process, pathname).insert(0,process.initSubset)
#         #if (options.mcversion == 'Summer12W0Jets'):
#         if (options.generator == 'w0jets'):
#             getattr(process, pathname).insert(0,process.filterWNJets)
#
#         getattr(process, pathname).insert(0,process.ak4PFCHSJetsBtag)
#         getattr(process, pathname).remove(process.goodOfflinePrimaryVertices)
#         getattr(process, pathname).insert(0,process.goodOfflinePrimaryVertices)
#         ## move the trigger to the beginning of the sequence
#         getattr(process, pathname).remove(process.hltFilter)
#         getattr(process, pathname).insert(0,process.hltFilter)
#
#     ## adjust lepton pT cut
#     import re
#     ptval="33" # adjust cut value here! [GeV]
#     relevantLeptonCollections = [process.tightElectronsEJ, process.goldenMuons]
#     exp = re.compile('(?:t\s?>\s?30)')
#     for lep in relevantLeptonCollections:
#         if(exp.search(lep.cut.pythonValue())!=None):
#             tmpExp=exp.sub("t > "+ptval, str(lep.cut.pythonValue()))
#             lep.cut=tmpExp
#             if(tmpExp.find("\'")>-1):
#                 lep.cut=tmpExp.strip("'")
