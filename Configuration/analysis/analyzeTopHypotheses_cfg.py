

#change 'lepton' in options to 'muon' or 'electron'

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('mcversion', 'RunIISpring16DR', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC campaign or data") #MC campaign is 'RunIISpring16DR' atm
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

options.register('taggerName', 'pfCombinedInclusiveSecondaryVertexV2BJetTags', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Tagger name in CMSSW")
options.register('csvm', 0.8, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "CSVM working point")  #0.89 on 74
options.register('nbjets', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Minimum number of bjets")

options.register('brCorrection', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Do BR correction (MadGraph)") #Off??? was on after nataliia merge
options.register('bSFNewRecipe', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Use new b-tag SF recipe")

options.register('addTriggerMatch' , False , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if trigger objects are matched to signal muons' )

options.register('cut', 'allCut', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Cut before kinFit") #if 'allCut' will aply all Cuts before the final Analyzer

options.register('dataElectronID', 'cutBasedElectronID-Spring15-25ns-V1', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "ElectronID configuration for data")  
options.register('mcElectronID', 'cutBasedElectronID-Spring15-25ns-V1', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "ElectronID configuration for MonteCarlo") 
#options.register('MCreHLT', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "set on True if running on MC reHLT")
options.register('MCtype', 'MC', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "set on MCreHLT if running on MC reHLT")

options.register('triggerless', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "set on True if no Triggerfilter is wanted i.e. for PU studies")
options.register('showLHEscaleWeightTypes', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "set on True to print out available scale Weight Types out of the LHEEventProducer")

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


allCut = (options.cut=='allCut') 
jet4Cut = (options.cut=='jet4Cut') 

# nTupel Output
#outputNTupel = cms.string( 'analyzeTop_{0}_{1}_{2}_Test.root'.format( options.lepton , options.lJesFactor , options.cut ) )
outputNTupel = cms.string('analyzeTop.root')

if (options.mcversion == "data"):
	data = True
	options.register('runOnMC', False , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if run on MC or real data' )
else:
	data = False
	options.register('runOnMC', True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, 'decide, if run on MC or real data' )
## top mass measurement
process = cms.Process("topMass")

if(data):
	electronID=options.dataElectronID
else:
	electronID=options.mcElectronID


## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.categories.append('TtSemiLeptonicEvent')
process.MessageLogger.cerr.TtSemiLeptonicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)


from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x import *

if data:  
	if (options.lepton=='muon'):
    		inputFiles = [
			'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/18383F36-2E1A-E611-8C57-02163E014186.root' 
			]
	else:
		inputFiles = [
			'/store/data/Run2016C/SingleElectron/MINIAOD/PromptReco-v2/000/276/062/00000/3A9AE0CC-8240-E611-BA98-02163E0144D1.root'
			]
else:
    inputFiles = [
      '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/007EB1F4-C028-E611-8AD7-0CC47A78A4A6.root'
      #'/store/mc/RunIISpring16MiniAODv2/TT_TuneEE5C_13TeV-powheg-herwigpp/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/D0B3C051-7F1B-E611-98F2-001EC9ADEA41.root'
      ]
    if (options.MCtype=='MCreHLT'):
      inputFiles = [
	'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/0064B539-803A-E611-BDEA-002590D0B060.root'
	]


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
triggerSelectionDataElectron  = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v* '
triggerSelectionMCElectron = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v* ' 

triggerSelectionDataMuon  = 'HLT_IsoMu22_v*' 
triggerSelectionDataMuon_Tk  = 'HLT_IsoTkMu22_v*' 
triggerSelectionMCMuon = 'HLT_IsoMu22_v*' 
triggerSelectionMCMuon_Tk = 'HLT_IsoTkMu22_v*' 
#triggerSelectionDataMuon = cms.vstring(
        #"HLT_IsoMu22_v*",
        #"HLT_IsoTkMu20_v*"
        ##"HLT_IsoMu24_eta2p1_v*",
        ##"HLT_IsoMu24_v*"
#)
if(options.triggerless):
	triggerSelectionDataElectron  = 'HLT* '
	triggerSelectionMCElectron = 'HLT* ' 
	triggerSelectionDataMuon  = 'HLT*' 
	triggerSelectionMCMuon = 'HLT*' 



#primary Vertex Cut? https://twiki.cern.ch/twiki/bin/view/CMS/TTbarHbbRun2ReferenceAnalysis_76XTransition  

# Step 1
#muonCut       = ''
#signalMuonCut = ''
#electronCut       = ''
#signalElectronCut = ''

#muonVertexMaxDZ = 0.5

# Step 2

# Step 3
useElecEAIsoCorr = options.useElecEAIsoCorr
useCalibElec     = options.useCalibElec

#electronGsfCut  =     '' 
#electronCalibCut = electronGsfCut.replace( 'ecalDrivenMomentum.', '' )

#electronGsfVetoCut  =     '' 
#electronCalibVetoCut = electronGsfVetoCut.replace( 'ecalDrivenMomentum.', '' )

#dileptonElectronVetoCut = ''
#conversionRejectionCut = ''

electronVetoCut = electronGsfVetoCut
electronCut = electronGsfCut
#signalElectronCut = 'pt > 30 && abs(eta) <2.1'

# Step 4 
#only interresting for the electron case
#conversionRejectioncCut=''

# Step 5

#jetCut = ''
#veryTightJetCut = cms.string('') 
#tightJetCut     = cms.string('') 
#looseJetCut     = cms.string('') 
#veryLooseJetCut = cms.string('') 



# Step 6
bTagSrc = 'selectedJets'
#bTagCut = 'bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.890' #0.814'
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
globalTagData = '80X_dataRun2_Prompt_ICHEP16JEC_v0' # Global tag: 80X_dataRun2_Prompt_v8 (2016B), 80X_dataRun2_Prompt_v9 (2016C & D), 
#usePrivateSQlite=False #do not use external JECs (sqlite file)
globalTagMC   = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
#globalTagData = 'DEFAULT'
usePrivateSQlite=True #use external JECs (sqlite file)

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

if (options.lepton=='muon'):
 triggerSelection       = triggerSelectionDataMuon

 triggerObjectSelection = triggerObjectSelectionData
 if runOnMC:
   triggerSelection       = triggerSelectionMCMuon
   triggerObjectSelection = triggerObjectSelectionMC
else:
  triggerSelection       = triggerSelectionDataElectron
  if runOnMC:
    triggerSelection       = triggerSelectionMCElectron
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
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") #FIXME this would be the recommended one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")


process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if runOnMC:
  if globalTagMC == 'DEFAULT':
    process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:run2_mc' )
  else:
    process.GlobalTag.globaltag = globalTagMC
else:
  if globalTagData == 'DEFAULT':
    process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:run2_data' )
  else:
    process.GlobalTag.globaltag = globalTagData

#override JEC   
if usePrivateSQlite: #now false
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import CondCore.CondDB.CondDB_cfi
    import os
    if runOnMC:
        era="Spring16_25nsV6_MC"  
    else:
        era="Spring16_25nsV6_DATA"   #TODO does he find this?
    dBFile = os.path.expandvars("$CMSSW_BASE/src/TopMass/data/JEC/"+era+".db")
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
    
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *
#from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJets = updatedPatJets.clone(
				  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC") )
				  )



#redo jets and MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#remove HF candidates
process.noHFCands = cms.EDFilter("CandPtrSelector",
                                     src=cms.InputTag("packedPFCandidates"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
if runOnMiniAOD:
	if data:
   		 runMetCorAndUncFromMiniAOD(process,
		   	        	            isData= (data),
		                	     pfCandColl=cms.InputTag("noHFCands"),#comment for MET with HF #FIXME noHFCands not recommended
                          		     jecUncFile='TopMass/data/JEC/Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt', 
                            	  	 )

	else:
   		 runMetCorAndUncFromMiniAOD(process,
		   	        	            isData= (data),
		                	        #pfCandColl=cms.InputTag("noHFCands"),#comment for MET with HF #FIXME noHFCands not recommended
                                    pfCandColl=cms.InputTag("packedPFCandidates"),
                                    #reclusterJets=True, #True needed for NoHF, but does not work... 
                                    #recoMetFromPFCs=True, #needed for NoHF, but kills Btagging....
                          		     jecUncFile='TopMass/data/JEC/Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt', 
                                    #jetCollUnskimmed="slimmedJets", #is default already...
                                    #jetColl="selectedPatJets"
                            	  	 )	
else:
    runMETCorrectionsAndUncertainties(process,
                                      isData= (data),
                                      pfCandColl=cms.InputTag("noHFCands"),#comment for MET with HF
                                      jecUncFile='TopMass/data/JEC/Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt', 
                                      postfix="NoHF")

del process.slimmedMETs.t01Variation #brute the swarm


process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
## Define PFJetsCHS
process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)

#add mnodules to redo the b-tagging on-the-fly using unscheduled mode
#from PhysicsTools.PatAlgos.tools.jetTools import *
## b-tag discriminators
bTagDiscriminators = [
    #'pfCombinedInclusiveSecondaryVertexV2BJetTags'
    options.taggerName
]

jetCorrectionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

if (data):
    jetCorrectionLevels.append('L2L3Residual')  #FIXME looks plausible like https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC


from PhysicsTools.PatAlgos.tools.jetTools import *
## call addJetCollection to get get b tagging modules prepared,  FIXME but not used anymore
addJetCollection(
    process,
    labelName = 'Updated',
    jetSource = cms.InputTag('ak4PFJetsCHS'), #btw what is 'CHS' difference from AK4 to AK8? -> explanaition: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetAnalysis 
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', jetCorrectionLevels, 'None'),
    genJetCollection = cms.InputTag('slimmedGenJets'),#cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #getJetMCFlavour = runOnMC, #TODO does not work at the moment
    getJetMCFlavour = False,
    algo = 'AK',
    rParam = 0.4
)


#Debug
#process.patJets.addJetCorrFactors = cms.bool(False)

#enable GenJets and GenPartons for the Output Jets
if runOnMC: 
    process.patJetsUpdated.addGenJetMatch=cms.bool(True) 
    process.patJetsUpdated.addGenPartonMatch=cms.bool(True) 


from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
if runOnMiniAOD:
    adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))


process.patJetsUpdated.addBTagInfo = cms.bool(True)
if runOnMC == False:
   process.patJetCorrFactorsUpdated.levels.append('L2L3Residual')
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
#  else:
#    if runOnMC:
#      from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
#      inputFiles = filesRelValProdTTbarAODSIM
#    else:
#      from PhysicsTools.PatAlgos.patInputFiles_cff import filesSingleMuRECO # not available at CERN
#      inputFiles = filesSingleMuRECO

process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" ) 

process.source.fileNames = inputFiles
# maximum number of events
process.maxEvents.input = maxEvents 


#TODO
runMetCorAndUncFromMiniAOD

###
### Selection configuration
###

# Individual steps

# Step 0

from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults

process.triggerSelection = triggerResults.clone( triggerConditions = [ triggerSelection ] )

from HLTrigger.HLTfilters.hltHighLevel_cfi import *

if (options.lepton=='muon'):
  import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
  if (options.MCtype=='MCreHLT'):
    process.filter_any_explicit = hlt.hltHighLevel.clone(
      TriggerResultsTag=cms.InputTag("TriggerResults","","HLT2"),
      HLTPaths = [triggerSelectionMCMuon, triggerSelectionMCMuon_Tk],
      throw = False
    )
  else:
    process.filter_any_explicit = hlt.hltHighLevel.clone(
      HLTPaths = [triggerSelectionDataMuon, triggerSelectionDataMuon_Tk],
      throw = False
    )
    if runOnMC:
      process.filter_any_explicit = hlt.hltHighLevel.clone(
	HLTPaths = [triggerSelectionMCMuon, triggerSelectionMCMuon_Tk],
	throw = False
      )
  process.sStandAloneTrigger = cms.Sequence( process.filter_any_explicit
					  )
else:  
  process.sStandAloneTrigger = cms.Sequence( process.triggerSelection
					  )
  process.pStandAloneTrigger = cms.Path( process.sStandAloneTrigger )

from TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cfi import metFiltersMiniAOD
process.metFiltersMiniAOD = metFiltersMiniAOD

process.eventCleaningMiniAOD = cms.Sequence(process.metFiltersMiniAOD)
process.eventCleaningMiniAODData = cms.Sequence()
process.eventCleaningMiniAODMC = cms.Sequence()
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

process.pStandAloneEventCleaning = cms.Path( process.sStandAloneEventCleaning )

from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
process.goodOfflinePrimaryVertices = goodOfflinePrimaryVertices.clone( filter = True ) #entspricht der PV selection empfelung aus dem Twiki (stand 76x)
if runOnMiniAOD:
  process.goodOfflinePrimaryVertices.src = 'offlineSlimmedPrimaryVertices'
process.sStandAloneGoodVertex = cms.Sequence( process.goodOfflinePrimaryVertices
                                            )
process.pStandAloneGoodVertex = cms.Path( process.sStandAloneGoodVertex )


# Step 1

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import selectedMuons, preSignalMuons, signalMuons, standAloneSignalMuonFilter, selectedElectrons, preSignalElectrons, signalElectrons, standAloneSignalElectronFilter

if (options.lepton=='muon'):
 process.selectedMuons = selectedMuons.clone( cut = muonCut )
 if runOnMiniAOD:
  	 process.selectedMuons.src = 'slimmedMuons'
 process.preSignalMuons = preSignalMuons.clone( cut = signalMuonCut )
 process.signalMuons = signalMuons.clone( maxDZ = muonVertexMaxDZ )
 if runOnMiniAOD:
  	 process.signalMuons.vertexSource = 'offlineSlimmedPrimaryVertices'
 process.standAloneSignalLeptonFilter = standAloneSignalMuonFilter.clone( ) 
else:
 process.selectedElectrons = selectedElectrons.clone( cut = (electronCut.format( electronID )) )  
 if runOnMiniAOD:
  	 process.selectedElectrons.src = 'slimmedElectrons'
 process.preSignalElectrons = preSignalElectrons.clone( cut = signalElectronCut ) 
 process.signalElectrons = signalElectrons.clone()  
 process.standAloneSignalLeptonFilter = standAloneSignalElectronFilter.clone( ) 

process.sStandAloneSignalLepton = cms.Sequence( process.standAloneSignalLeptonFilter )
process.pStandAloneSignalLepton = cms.Path( process.sStandAloneSignalLepton )

# Step 2

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import standAloneLooseMuonVetoFilter, LooseMuonVetoFilter

if (options.lepton=='muon'):
 process.standAloneLooseMuonVetoFilter = standAloneLooseMuonVetoFilter.clone()
 process.sStandAloneLooseMuonVeto = cms.Sequence( process.standAloneLooseMuonVetoFilter )
else:
 process.selectedMuons = selectedMuons.clone( cut = muonCut )
 if runOnMiniAOD:
 	  process.selectedMuons.src = 'slimmedMuons'
 process.preSignalMuons = preSignalMuons.clone( cut = signalMuonCut )
 process.signalMuons = signalMuons.clone( maxDZ = muonVertexMaxDZ )
 if runOnMiniAOD:
 	  process.signalMuons.vertexSource = 'offlineSlimmedPrimaryVertices'
 process.LooseMuonVetoFilter = LooseMuonVetoFilter.clone()
 process.sStandAloneLooseMuonVeto = cms.Sequence( process.LooseMuonVetoFilter )

process.pStandAloneLooseMuonVeto = cms.Path( process.sStandAloneLooseMuonVeto )


# Step 3

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
  electronGsfVetoCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )
  electronCalibVetoCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )
  electronVetoCut.replace( '-1.0*userIsolation("User1Iso")', '-0.5*puChargedHadronIso' )

if useCalibElec: #not used atm (option at the document beginning)
  from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import electronsWithRegression, calibratedElectrons
  process.electronsWithRegression = electronsWithRegression.clone()
  if runOnMiniAOD:
    if useElecEAIsoCorr:
      process.electronsWithRegression.inputElectronsTag = 'electronsWithEA03Iso'
    else:
      process.electronsWithRegression.inputElectronsTag = 'slimmedElectrons'
    process.electronsWithRegression.vertexCollection  = 'offlineSlimmedPrimaryVertices'
  process.calibratedElectrons = calibratedElectrons.clone( isMC = runOnMC )
  if runOnMC:
    process.calibratedElectrons.inputDataset = 'Summer12_LegacyPaper' #not used atm
  else:
    process.calibratedElectrons.inputDataset = '22Jan2013ReReco' #not used atm
  process.RandomNumberGeneratorService = cms.Service( "RandomNumberGeneratorService"
                                                    , calibratedElectrons = cms.PSet( initialSeed = cms.untracked.uint32( 1 )
                                                                                    , engineName  = cms.untracked.string('TRandom3')
                                                                                    )
                                                    )
  electronVetoCut = electronCalibVetoCut



from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import selectedElectrons, standAloneElectronVetoFilter, dileptonElectronVetoFilter, selectedElectrons4Veto
if (options.lepton=='muon'):
 process.selectedElectrons = selectedElectrons.clone( cut = (electronVetoCut.format( electronID )) )
 if useCalibElec:
 	 process.selectedElectrons.src = 'calibratedElectrons'
 elif useElecEAIsoCorr and runOnMiniAOD:
  	 process.selectedElectrons.src = 'electronsWithEA03Iso'
 elif runOnMiniAOD:
 	 process.selectedElectrons.src = 'slimmedElectrons'
 process.ElectronVetoFilter = standAloneElectronVetoFilter.clone()

else:
 process.selectedElectrons4Veto = selectedElectrons4Veto.clone( cut = (dileptonElectronVetoCut.format( electronID )) )
 if useCalibElec:
  	 process.selectedElectrons4Veto.src = 'calibratedElectrons'
 elif useElecEAIsoCorr and runOnMiniAOD:
  	 process.selectedElectrons4Veto.src = 'electronsWithEA03Iso'
 elif runOnMiniAOD:
  	 process.selectedElectrons4Veto.src = 'slimmedElectrons'
 process.ElectronVetoFilter = dileptonElectronVetoFilter.clone()
 
process.sElectronVetoX = cms.Sequence( process.ElectronVetoFilter )
process.pElectronVetoX = cms.Path( process.sElectronVetoX )

# Step 4
if (options.lepton=='electron'):
 from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import conversionRejectionFilter
 process.conversionRejectionFilter = conversionRejectionFilter.clone( cut = (conversionRejectionCut.format( electronID )) )
 process.sConversionRejectionFilter = cms.Sequence( process.conversionRejectionFilter )
 process.pConversionRejectionFilter = cms.Path( process.sConversionRejectionFilter )

# Step 5

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import selectedJets
process.selectedJets = selectedJets.clone( cut = jetCut )

process.goodOfflinePrimaryVertices.taggedMode=cms.untracked.bool( True )

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import signalVeryTightJets, standAloneSignalVeryTightJetsFilter
process.signalVeryTightJets = signalVeryTightJets.clone( cut = veryTightJetCut )
process.standAloneSignalVeryTightJetsFilter = standAloneSignalVeryTightJetsFilter.clone(  ) 
process.sStandAlone1Jet = cms.Sequence( process.standAloneSignalVeryTightJetsFilter )
process.pStandAlone1Jet = cms.Path( process.sStandAlone1Jet )

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import signalTightJets, standAloneSignalTightJetsFilter
process.signalTightJets = signalTightJets.clone( cut = tightJetCut )
process.standAloneSignalTightJetsFilter = standAloneSignalTightJetsFilter.clone( )
process.sStandAlone2Jets = cms.Sequence( process.standAloneSignalTightJetsFilter )
process.pStandAlone2Jets = cms.Path( process.sStandAlone2Jets )

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import signalLooseJets, standAloneSignalLooseJetsFilter
process.signalLooseJets = signalLooseJets.clone( cut = looseJetCut )
process.standAloneSignalLooseJetsFilter = standAloneSignalLooseJetsFilter.clone(  )
process.sStandAlone3Jets = cms.Sequence( process.standAloneSignalLooseJetsFilter )
process.pStandAlone3Jets = cms.Path( process.sStandAlone3Jets )

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import signalVeryLooseJets, standAloneSignalVeryLooseJetsFilter
process.signalVeryLooseJets = signalVeryLooseJets.clone( cut = veryLooseJetCut )
process.standAloneSignalVeryLooseJetsFilter = standAloneSignalVeryLooseJetsFilter.clone(  )
process.sStandAlone4Jets = cms.Sequence( process.standAloneSignalVeryLooseJetsFilter )
process.pStandAlone4Jets = cms.Path( process.sStandAlone4Jets )

if (options.lepton=='electron'):
    process.cleanedJets = cms.EDProducer("PATJetCleaner", 
                                        src = cms.InputTag("selectedJets"),
                                         # preselection (any string-based cut on pat::Jet)     
                                         preselection = cms.string(''), 
                                         # overlap checking configurables
                                         checkOverlaps = cms.PSet(
                                                                  electrons = 
                                                                  cms.PSet(
                                                                       src       = cms.InputTag("signalElectrons"), 
                                                                       algorithm = cms.string("byDeltaR"),
                                                                       preselection        = cms.string(""),
                                                                       deltaR              = cms.double(0.3),
                                                                       checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
                                                                       pairCut             = cms.string(""),
                                                                       requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
                                                                       )
                                                              ),
                                         # finalCut (any string-based cut on pat::Jet)
                                         finalCut = cms.string(''))


else:
    process.cleanedJets = cms.EDProducer("PATJetCleaner", 
                                         src = cms.InputTag("selectedJets"),
                                         # preselection (any string-based cut on pat::Jet)     
                                         preselection = cms.string(''), 
                                         # overlap checking configurables
                                         checkOverlaps = cms.PSet(muons = 
                                                                cms.PSet(
                                                                         src       = cms.InputTag("signalMuons"),
                                                                         algorithm = cms.string("byDeltaR"),
                                                                         preselection        = cms.string(""), 
                                                                         deltaR              = cms.double(0.3),
                                                                         checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref          
                                                                         pairCut             = cms.string(""),
                                                                         requireNoOverlaps   = cms.bool(True)
                                                                         )
                                                              ),
                                         # finalCut (any string-based cut on pat::Jet)
                                         finalCut = cms.string(''))

process.signalVeryLooseJets.src = 'cleanedJets'


# Step 6

from TopMass.Configuration.patRefSel_refElectronJets_refMuJets_76x_cfi import selectedBTagJets, standAloneSignalBTagsFilter
process.selectedBTagJets = selectedBTagJets.clone( src = bTagSrc
                                                 , cut = bTagCut
                                                 )
process.standAloneSignalBTagsFilter = standAloneSignalBTagsFilter.clone( minNumber = minBTags)
process.sStandAloneBTags = cms.Sequence( process.standAloneSignalBTagsFilter )
process.pStandAloneBTags = cms.Path( process.sStandAloneBTags )



#DummyFilter
#A Producer to read out, which EDFilters were passed
from TopMass.TopEventTree.FilterDummy_cfi import FilterDummy
process.FilterDummyTrigger = FilterDummy.clone()
process.FilterDummyEventCleaning = FilterDummy.clone()
process.FilterDummyGoodVertex = FilterDummy.clone()
process.FilterDummySignalLepton = FilterDummy.clone( )
process.FilterDummyLooseMuonVeto = FilterDummy.clone()
process.FilterDummyElectronVeto = FilterDummy.clone()
process.FilterDummyConversionRejection = FilterDummy.clone()
process.FilterDummy1Jet = FilterDummy.clone( )
process.FilterDummy2Jets = FilterDummy.clone()
process.FilterDummy3Jets = FilterDummy.clone()
process.FilterDummy4Jets = FilterDummy.clone()
process.FilterDummyBTags = FilterDummy.clone()

#for Debugging
#from TopMass.TopEventTree.JetDebug_cfi import JetDebug
#process.JD1 = JetDebug.clone(jets = cms.InputTag("patJets"), evtSolLabel= cms.InputTag("ttSemiLepEvent") )


# Trigger matching

if addTriggerMatch:
 if (options.lepton=='muon'):  #FIXME electron case missing, bot available in TopQuarkAnalysis.Configuration.patRefSel_triggerMatching ...
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
    #                               , 'keep *_signalMuonsTriggerMatch_*_*'
    #                               ]
    # process.out.outputCommands += [ 'drop *_signalMuons_*_*'
    #                               , 'keep *_signalMuonsTriggerMatch_*_*'
    #                               ]






## sequence for ttGenEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")


## sequence for TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")

## enable additional per-event printout from the TtSemiLeptonicEvent
process.ttSemiLepEvent.verbosity = 0

#new 7.3.16
process.decaySubset.runMode = cms.string("Run2")

## TRIGGER
from HLTrigger.HLTfilters.hltHighLevel_cfi import *

if os.getenv('CMSSW_VERSION').startswith('CMSSW_8_1_'):  #TODO do we want this on 76?
  if (options.MCtype=='MCreHLT'):
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT2", throw=True)
    process.triggerSelection.hltResults = cms.InputTag("TriggerResults","","HLT2")
  else:
    process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", throw=True)
    process.triggerSelection.hltResults = cms.InputTag("TriggerResults","","HLT")
    
## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
if (options.lepton=='muon'):
 setForAllTtSemiLepHypotheses(process, "leps", "signalMuons")
else:
 if (options.mcversion=='genLevel'):
  setForAllTtSemiLepHypotheses(process, "leps", 'goodElectronsEJ')
 else:
  useElectronsForAllTtSemiLepHypotheses(process, 'signalElectrons')
setForAllTtSemiLepHypotheses(process, "jets", "signalVeryLooseJets")
setForAllTtSemiLepHypotheses(process, "maxNJets", 4)
if (options.uncFactor!=1.0) :
 setForAllTtSemiLepHypotheses(process, "mets", "scaledMET:scaledMETs")
setForAllTtSemiLepHypotheses(process, "mets", "patPFMetT1")


setForAllTtSemiLepHypotheses(process, "maxNComb", -1)
if data:
    #setForAllTtSemiLepHypotheses(process, "mets", "patMETs")
    setForAllTtSemiLepHypotheses(process, "jetCorrectionLevel", "L2L3Residual")  #FIXME what does this do? was commented out before gives Exception Message: This JEC level L2L3Residual does not exist.


## change jet-parton matching algorithm #FIXME enstspricht nicht der empfehlung in seinem Kommentar
process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
process.ttSemiLepJetPartonMatch.maxNJets  = -1
process.ttSemiLepJetPartonMatch.useMaxDist = True

# consider b-tagging in event reconstruction
process.hitFitTtSemiLepEventHypothesis.bTagAlgo = options.taggerName
process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = options.csvm
process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = options.csvm
#process.hitFitTtSemiLepEventHypothesis.minBDiscBJets     = 1.0
#process.hitFitTtSemiLepEventHypothesis.maxBDiscLightJets = 3.0
process.hitFitTtSemiLepEventHypothesis.useBTagging       = True

addTtSemiLepHypotheses(process,
                       ["kHitFit", "kMVADisc"]
                       )

if (data): removeTtSemiLepHypGenMatch(process)

## load HypothesisAnalyzer
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import analyzeHypothesis
process.analyzeHitFit = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key", jets= "signalVeryLooseJets")
if (options.lepton=="electron"): 
     process.analyzeHitFit.lepton = 11
from TopMass.TopEventTree.JetEventAnalyzer_cfi import analyzeJets
process.analyzeJets = analyzeJets.clone(jets = "signalVeryLooseJets",met = "patPFMetT1")
from TopMass.TopEventTree.WeightEventAnalyzer_cfi import analyzeWeights
process.analyzeWeights = analyzeWeights.clone(

                                              mcWeight        = options.mcWeight,
					                          puSrc     = cms.InputTag("slimmedAddPileupInfo"),
                                              puWeightSrc     = cms.InputTag("eventWeightPUsysNo"  , "eventWeightPU"),
                                              puWeightUpSrc   = cms.InputTag("eventWeightPUsysUp"  , "eventWeightPU"),
                                              puWeightDownSrc = cms.InputTag("eventWeightPUsysDown", "eventWeightPU"),
                                              savePDFWeights = True,
                                              brCorrection   = options.brCorrection
                                             )
if(options.showLHEscaleWeightTypes):
    process.analyzeWeights.showLHEweightTypes=True

if (options.lepton=='electron'):
    
	process.analyzeWeights.lepIDWeightSrc      = cms.InputTag("effSFElectronEventWeightLepID")
	process.analyzeWeights.lepIDWeightSrcUp    = cms.InputTag("effSFElectronEventWeightUpLepID")
	process.analyzeWeights.lepIDWeightSrcDown  = cms.InputTag("effSFElectronEventWeightDownLepID")
	
	process.analyzeWeights.isoWeightSrc      = cms.InputTag("effSFElectronEventWeightIso")
	process.analyzeWeights.isoWeightSrcUp    = cms.InputTag("effSFElectronEventWeightUpIso")
	process.analyzeWeights.isoWeightSrcDown  = cms.InputTag("effSFElectronEventWeightDownIso")
	
	process.analyzeWeights.triggerWeightSrc      = cms.InputTag("effSFElectronEventWeightTrigger")
	process.analyzeWeights.triggerWeightSrcUp    = cms.InputTag("effSFElectronEventWeightUpTrigger")
	process.analyzeWeights.triggerWeightSrcDown  = cms.InputTag("effSFElectronEventWeightDownTrigger")
	

if (options.lepton=='muon'):

	process.analyzeWeights.lepIDWeightSrc      = cms.InputTag("effSFMuonEventWeightLepID")
	process.analyzeWeights.lepIDWeightSrcUp    = cms.InputTag("effSFMuonEventWeightUpLepID")
	process.analyzeWeights.lepIDWeightSrcDown  = cms.InputTag("effSFMuonEventWeightDownLepID")
	
	process.analyzeWeights.isoWeightSrc      = cms.InputTag("effSFMuonEventWeightIso")
	process.analyzeWeights.isoWeightSrcUp    = cms.InputTag("effSFMuonEventWeightUpIso")
	process.analyzeWeights.isoWeightSrcDown  = cms.InputTag("effSFMuonEventWeightDownIso")
	
	process.analyzeWeights.triggerWeightSrc      = cms.InputTag("effSFMuonEventWeightTrigger")
	process.analyzeWeights.triggerWeightSrcUp    = cms.InputTag("effSFMuonEventWeightUpTrigger")
	process.analyzeWeights.triggerWeightSrcDown  = cms.InputTag("effSFMuonEventWeightDownTrigger")
	
process.ttSemiLepHypGenMatch.useBReg = cms.bool(False)
process.ttSemiLepHypMVADisc.useBReg = cms.bool(False)

if runOnMiniAOD:
	process.initSubset.src = "prunedGenParticles"
	process.decaySubset.src = "prunedGenParticles"

#process.analyzeHitFitSel = process.analyzeHitFit.clone(topBranchName = 'topSel.')
#process.analyzer         = cms.Sequence(process.makeGenEvt+process.makeTtSemiLepEvent+process.analyzeHitFit+process.analyzeJets+process.analyzeWeights )
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.analyzer         = cms.Sequence(process.analyzeHitFit+process.analyzeJets+process.analyzeWeights) 

if runOnMC:
    #use smeared jets and MET
    process.analyzeJets.met = "patPFMetT1Smear"
    setForAllTtSemiLepHypotheses(process, "mets", "patPFMetT1Smear")
    #process.cleanedJets.src = src = "patSmearedJets" #why was that?
    
#atm first scaling, then cleaning!
   #*****************************************************************************
   #add systematic variations here using products from the MET uncertainty tool
   #*****************************************************************************
    ## configure JetEnergyScale tool 
process.load("TopAnalysis.TopUtils.JetEnergyScale_cff")
from TopAnalysis.TopUtils.JetEnergyScale_cff import *

scaledJetEnergy.scaleType    = options.scaleType
scaledJetEnergy.JECUncSrcFile= "TopMass/data/JEC/Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt" #Uncertainty sources 
scaledJetEnergy.sourceName   = options.jessource
scaledJetEnergy.flavor       = options.flavor
scaledJetEnergy.scaleFactor  = options.lJesFactor
scaledJetEnergy.scaleFactorB = options.bJesFactor

   #see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Uncertainty
if runOnMiniAOD:
    #resolutionNominal = [1.095, 1.120, 1.097, 1.103, 1.118, 1.100, 1.162, 1.160, 1.161, 1.209, 1.564, 1.384, 1.216] #76X
    #resolutionNominal = [1.122, 1.167, 1.168, 1.029, 1.115, 1.041, 1.167, 1.094, 1.168, 1.266, 1.595, 0.998, 1.226]
    resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
elif options.mcversion == "genLevel":
    resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
else:
    resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#elif options.mcversion == "genLevel":
#    resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#else:
#    resolutionNominal = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   #[1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056] #FIXME alle auf 1.0?  nmbrs from Run I
#resolutionUnc     = [0.018, 0.028,0.017,0.033,0.014,0.033,0.044,0.048,0.060,0.059,0.321,0.033,0.050]#76X
resolutionUnc     = [0.026, 0.048,0.046,0.066,0.03,0.062,0.086,0.093,0.120,0.132,0.175,0.066,0.145]
	
scaledJetEnergy.resolutionFactors   = [] 
for nom,unc in zip(resolutionNominal,resolutionUnc):
  scaledJetEnergy.resolutionFactors.append(nom+options.resolutionSigma*unc)
	
#scaledJetEnergy.resolutionEtaRanges   = [0.0,0.5, 0.5,1.1, 1.1,1.7, 1.7,2.3, 2.3,2.8, 2.8,3.2, 3.2,-1.]
scaledJetEnergy.resolutionEtaRanges   = [0.0,0.5, 0.5,0.8, 0.8,1.1, 1.1,1.3, 1.3,1.7, 1.7,1.9, 1.9,2.1, 2.1,2.3, 2.3,2.5, 2.5,2.8, 2.8,3.0, 3.0,3.2, 3.2,-1.]

scaledJetEnergy.inputJets    = "patJets" #"patJetsUpdated" 
scaledJetEnergy.inputMETs    = "patPFMetT1" 

process.selectedJets.src ="scaledJetEnergy:patJets"  #"scaledJetEnergy:patJetsUpdated" 
#soll das atm verwendet werden? uncertainties from http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-15-003/index.html
if (options.lepton=='electron' and options.eesShift!=0.0):
	## electron shift
	process.load("TopAnalysis.TopUtils.ElectronEnergyScale_cfi")
	from TopAnalysis.TopUtils.ElectronEnergyScale_cfi import *
	
	process.scaledElectronEnergy.src      = "selectedElectrons" 
	process.scaledElectronEnergy.mets     = "scaledJetEnergy:patPFMetT1"
   # process.scaledElectronEnergy.binning  = cms.VPSet(
   #                   cms.PSet(
   #                       binSelection = cms.string('isEB'),
   #                       binUncertainty = cms.double(0.01)
   #                   ),
   #                   cms.PSet(
   #                       binSelection = cms.string('!isEB'),
   #                       binUncertainty = cms.double(0.01)
   #                   )
   #                                                   )
	process.scaledElectronEnergy.shiftBy  = options.eesShift
	process.preSignalElectrons.src   = "scaledElectronEnergy:selectedElectrons" 
elif(options.mesShift!=0.0):
	## muon shift
	process.load("TopAnalysis.TopUtils.MuonEnergyScale_cfi")
	from TopAnalysis.TopUtils.MuonEnergyScale_cfi import *
	
	process.scaledMuonEnergy.src      = "selectedMuons"  
	process.scaledMuonEnergy.mets     = "scaledJetEnergy:patPFMetT1" 
        process.scaledMuonEnergy.uncertainty = cms.double(0.005)
	process.scaledMuonEnergy.shiftBy  = options.mesShift
	process.preSignalMuons.src   = "scaledMuonEnergy:selectedMuons" 

## unclustered energy scale
process.load("TopAnalysis.TopUtils.UnclusteredMETScale_cfi")
from TopAnalysis.TopUtils.UnclusteredMETScale_cfi import *

process.scaledMET.inputJets       = "scaledJetEnergy:selectedPatJets"   #FIXME selectedPatJetsUpdated ???
process.scaledMET.inputMETs       = "scaledMuonEnergy:METs"
process.scaledMET.inputElectrons  = "scaledElectronEnergy:selectedElectrons"
process.scaledMET.inputMuons      = "scaledMuonEnergy:selectedPatMuons"
process.scaledMET.scaleFactor     = options.uncFactor



process.makeJets = cms.Path(process.patJets)






 
#process.RAWSIMoutput = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('output_Test.root') )
#process.epath = cms.EndPath(process.RAWSIMoutput)

##Good for debugging
#process.evtDump = cms.EDAnalyzer("DumpEvent")
#process.dumpper = cms.Sequence(process.evtDump)


process.sTrigger       = cms.Sequence( process.sStandAloneTrigger + process.FilterDummyTrigger 
                                     )
process.sEventCleaning = cms.Sequence( process.sStandAloneEventCleaning + process.FilterDummyEventCleaning
                                     )
process.sGoodVertex    = cms.Sequence( process.sStandAloneGoodVertex + process.FilterDummyGoodVertex
                                     )
process.sSignalLepton    = cms.Sequence( process.sStandAloneSignalLepton + process.FilterDummySignalLepton
                                     )
process.sLooseMuonVeto = cms.Sequence( process.sStandAloneLooseMuonVeto + process.FilterDummyLooseMuonVeto
                                     )
process.sElectronVeto  = cms.Sequence( process.sElectronVetoX + process.FilterDummyElectronVeto
                                     )
if (options.lepton=='electron'):
 process.sConversionRejection = cms.Sequence( process.sConversionRejectionFilter + process.FilterDummyConversionRejection 
				     )
else:
 process.sConversionRejection = cms.Sequence( process.FilterDummyConversionRejection 
				     )
process.s1Jet          = cms.Sequence( process.sStandAlone1Jet + process.FilterDummy1Jet
                                     )
process.s2Jets         = cms.Sequence( process.sStandAlone2Jets + process.FilterDummy2Jets
                                     )
process.s3Jets         = cms.Sequence( process.sStandAlone3Jets + process.FilterDummy3Jets
                                     )
process.s4Jets         = cms.Sequence( process.sStandAlone4Jets + process.FilterDummy4Jets
                                     )
process.sBTags         = cms.Sequence( process.sStandAloneBTags + process.FilterDummyBTags
                                     )

process.serialFilter = cms.Sequence()



if allCut:
	process.serialFilter = cms.Sequence( #process.sTrigger
					process.sEventCleaning
					+process.sGoodVertex 
					+process.sSignalLepton
					+process.sLooseMuonVeto
					+process.sElectronVeto
					+process.sConversionRejection
					+process.s1Jet
					+process.s2Jets
					+process.s3Jets
					+process.s4Jets
					+process.sBTags
					)
else:
	if jet4Cut:
		process.serialFilter = cms.Sequence( #process.sTrigger
					process.sEventCleaning
					+process.sGoodVertex 
					+process.sSignalLepton
					+process.sLooseMuonVeto
					+process.sElectronVeto
					+process.sConversionRejection
					+process.s1Jet
					+process.s2Jets
					+process.s3Jets
					+process.s4Jets
					)		
	#else:
	#	process.serialFilter = cms.Sequence(process.sTrigger)

process.pTrigger       = cms.Path( process.sTrigger )
process.pEventCleaning = cms.Path( process.sEventCleaning )
process.pGoodVertex    = cms.Path( process.sGoodVertex)
process.pSignalLepton    = cms.Path( process.sSignalLepton)
process.pLooseMuonVeto = cms.Path( process.sLooseMuonVeto )
process.pElectronVeto  = cms.Path( process.sElectronVeto )
process.pConversionRejection = cms.Path (process.sConversionRejection )
process.p1Jet          = cms.Path( process.s1Jet )
process.p2Jets         = cms.Path( process.s2Jets)
process.p3Jets         = cms.Path( process.s3Jets)
process.p4Jets         = cms.Path( process.s4Jets )
process.pBTags         = cms.Path( process.sBTags )


from TopMass.TopEventTree.FilterDummyAnalyzer_cfi import FilterDummyAnalyzer
process.FDAnalyzer = FilterDummyAnalyzer.clone()


if not data:
  # calculate b tag efficiencies
  process.load("TopAnalysis.TopAnalyzer.BTagEfficiencyAnalyzer_cfi")
  # NOTE: process needs to be named bTagEff, so that BTagSFEventWeight.cc can find the histo
  process.bTagEff = process.analyzeBTagEfficiency.clone(jets         = "signalVeryLooseJets",
							bTagAlgo    = options.taggerName,
							bTagDiscrCut = options.csvm, ## CSVM
							weight       = process.analyzeWeights.puWeightSrc
						      )
  process.pbTagEff = cms.Path(process.bTagEff)
  
process.pAnalysis = cms.Path(process.serialFilter + process.analyzer + process.FDAnalyzer) 

  #process.schedule = cms.Schedule(
    #process.makeJets,   
    #process.pTrigger,
    #process.pEventCleaning,
    #process.pGoodVertex,
    #process.pSignalLepton,
    #process.pLooseMuonVeto,
    #process.pElectronVeto,
    #process.pConversionRejection,
    #process.p1Jet,
    #process.p2Jets,
    #process.p3Jets,
    #process.p4Jets,
    #process.pbTagEff,
    #process.pBTags,
    #process.pAnalysis
    #) 
 


process.schedule = cms.Schedule(
  process.makeJets,   
  process.pTrigger,
  process.pEventCleaning,
  process.pGoodVertex,
  process.pSignalLepton,
  process.pLooseMuonVeto,
  process.pElectronVeto,
  process.pConversionRejection,
  process.p1Jet,
  process.p2Jets,
  process.p3Jets,
  process.p4Jets,
  process.pbTagEff,
  process.pBTags,
  process.pAnalysis
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

        process.eventWeightPU.MCSampleFile             = "TopMass/data/PUhists/Run2_MC_2016B_asymptotic-v2_truePU.root"
        #process.eventWeightPU.MCSampleFile             = "TopMass/data/PUhists/Run2_MC_2016B_asymptotic-v2_truePU_reHLT.root"

    elif (options.mcversion == "RunIISpring16DR"):
	process.eventWeightPU.PUSource                 = cms.InputTag("slimmedAddPileupInfo")
        process.eventWeightPU.MCSampleHistoName        = "puhisto"
        process.eventWeightPU.DataHistoName            = "pileup"
        #process.eventWeightPU.MCSampleFile             = "TopMass/data/PUhists/Run2_MC_2016B_asymptotic-v2_truePU_reHLT.root"   #works for asymptotic_v12 PU generators
        process.eventWeightPU.MCSampleFile             = "TopMass/data/PUhists/Run2_MC_2016B_asymptotic-v2_truePU.root"   #works for asymptotic_v12 PU generators    
	process.eventWeightPU.DataFile		       = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS69mb.root" # --minBiasXsec 69mb;Mean 42.272160, RMS 10.574250
        process.eventWeightPUUp = process.eventWeightPU.clone(DataFile  = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS80mb.root") # --minBiasXsec 80mb; Mean 49.011200, RMS 12.260000
        process.eventWeightPUDown = process.eventWeightPU.clone(DataFile = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS58mb.root")# --minBiasXsec 58mb; Mean 49.011200, RMS 12.260000
    if (options.generator == "rd"):
        process.eventWeightPU.MCSampleHistoName        = ""
        process.eventWeightPU.DataHistoName            = ""
    

    process.eventWeightPUsysNo   = process.eventWeightPU.clone()
    process.eventWeightPUsysUp   = process.eventWeightPUUp.clone()
    process.eventWeightPUsysDown = process.eventWeightPUDown.clone()

    #### Parameters 'CreateWeight3DHisto' and 'Weight3DHistoFile' required for cff-file, but actually not used for Fall11 samples

    #### Configuration for Nominal PU Weights #TODO aktualisieren sollte passen

    process.eventWeightPUsysNo.WeightName          = "eventWeightPU"
    process.eventWeightPUsysNo.DataFile            = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS69mb.root"
    process.eventWeightPUsysNo.CreateWeight3DHisto = False
    #process.eventWeightPUsysNo.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3D.root"

    #### Configuration for PU Up Variations
    if(options.mcversion != "RunIISpring16DR"):
    	 process.eventWeightPUsysUp.WeightName          = "eventWeightPUUp"
    	 process.eventWeightPUsysUp.DataFile            = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS69mb.root"  
    	 process.eventWeightPUsysUp.CreateWeight3DHisto = False
    #process.eventWeightPUsysUp.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3DUp.root"

    #### Configuration for PU Down Variations

   	 process.eventWeightPUsysDown.WeightName          = "eventWeightPUDown"
    	 process.eventWeightPUsysDown.DataFile            = "TopMass/data/PUhists/Data_PUHist_Run2016B_12918pb_XS69mb.root"  
    	 process.eventWeightPUsysDown.CreateWeight3DHisto = False
    #process.eventWeightPUsysDown.Weight3DHistoFile   = "TopMass/Configuration/data/Weight3DDown.root"

    #### event weight sequence

    	 process.makeEventWeightsPU = cms.Sequence(process.eventWeightPUsysNo ) # * #atm no Variations
                                             # process.eventWeightPUsysUp   *
                                              #process.eventWeightPUsysDown  )

    ## ---
    ##    MC B-tag reweighting   
    ## ---
    ## load BTV database
    #process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
    #process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
    from RecoBTag.PerformanceDB.measure.Btag_mistag110711 import *
    from RecoBTag.PerformanceDB.measure.Btag_btag110711 import *
    process.load("TopAnalysis.TopUtils.BTagSFEventWeight_cfi")
    process.bTagSFEventWeight.jets     = "signalTightJets"
    process.bTagSFEventWeight.csv_filename = "TopMass/data/BTagSFs/CSVv2_ichep.csv"
    process.bTagSFEventWeight.bTagAlgo = options.taggerName #for CSVv2	
    process.bTagSFEventWeight.discrCut = options.csvm # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
    process.bTagSFEventWeight.version  = "CSVv2" # JP; CSVv2; cMVAv2 possible;
    if options.bSFNewRecipe:
        process.bTagSFEventWeight.newRecipe= True
        process.bTagSFEventWeight.maxJets  = 4
    process.bTagSFEventWeight.sysVar   = "" #bTagSFUp, bTagSFDown, bTagCjetSFUp, bTagCjetSFDown, misTagSFUp, misTagSFDown possible;
    process.bTagSFEventWeight.filename = "TopMass/data/BTagSFs/eff_bTag_2016B.root"  # NOTE: efficiency file should be up-to-date for each data sample 
    process.bTagSFEventWeight.verbose  = 0

    process.bTagSFEventWeightBTagSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
    process.bTagSFEventWeightBTagSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
    process.bTagSFEventWeightBTagCjetSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagCjetSFUp")
    process.bTagSFEventWeightBTagCjetSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagCjetSFDown")
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
        fragTargetFile = "TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star_rbLEP.root"  #TODO fragTargetSource to rbLEP
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
    ## scale factor for trigger and lepton selection efficiency  #TODO where is this used? as trigger weight in the WeightAnalyzer used as TriggerWeight, Lepton Weights label names imported but not used
    process.load("TopAnalysis.TopUtils.EffSFLepton2DEventWeight_cfi")
    process.effSFMuonEventWeight=process.effSFLepton2DEventWeight.clone()
    process.effSFMuonEventWeight.particles=cms.InputTag("signalMuons")
    process.effSFMuonEventWeight.sysVar   = cms.string("") #"noSys" ?
    #process.effSFMuonEventWeight.filename= "TopAnalysis/Configuration/data/MuonEffSF2D2012.root" #new ones from twiki 76x  page
    
    process.effSFMuonEventWeightLepID     = process.effSFMuonEventWeight.clone(filename="TopMass/data/LeptonSFs/MuonID_Z_RunBCD_prompt80X_7p65.root")
    process.effSFMuonEventWeightLepID.histoname="MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/pt_abseta_ratio"  #pt>120 is missing, solved via WeightEventAnalyzer
    
    process.effSFMuonEventWeightUpLepID  =process.effSFMuonEventWeightLepID.clone(sysVar = "combinedEffSFNormUpStat")
    process.effSFMuonEventWeightDownLepID=process.effSFMuonEventWeightLepID.clone(sysVar = "combinedEffSFNormDownStat")
    
    process.effSFMuonEventWeightIso     = process.effSFMuonEventWeight.clone(filename="TopMass/data/LeptonSFs/MuonIso_Z_RunBCD_prompt80X_7p65.root")
    process.effSFMuonEventWeightIso.histoname="MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/pt_abseta_ratio"  #pt>120 is missing, solved via WeightEventAnalyzer
    #process.effSFMuonEventWeightIso.verbose=cms.int32(15) #debugging 
    
    process.effSFMuonEventWeightUpIso  =process.effSFMuonEventWeightIso.clone(sysVar = "combinedEffSFNormUpStat")
    process.effSFMuonEventWeightDownIso=process.effSFMuonEventWeightIso.clone(sysVar = "combinedEffSFNormDownStat")
    
    
    process.effSFMuonEventWeightTrigger    = process.effSFMuonEventWeight.clone(filename="TopMass/data/LeptonSFs/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root")
    #process.effSFMuonEventWeight.filename="TopMass/data/LeptonSFs/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root"  #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHbbRun2ReferenceAnalysis_76XTransition
    process.effSFMuonEventWeightTrigger.histoname="IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097/efficienciesDATA/pt_abseta_DATA"  #pt>120 is missing, solved via WeightEventAnalyzer
    process.effSFMuonEventWeightTrigger.additionalSystErr=0.005  #+ 0.5% for single muon triggers to be applyied on top of statistical errors from root/pkl files
    process.effSFMuonEventWeightTrigger.verbose=cms.int32(0) 
    
    process.effSFMuonEventWeightUpTrigger  =process.effSFMuonEventWeightTrigger.clone(sysVar = "combinedEffSFNormUpStat")
    process.effSFMuonEventWeightDownTrigger=process.effSFMuonEventWeightTrigger.clone(sysVar = "combinedEffSFNormDownStat")





    process.effSFElectronEventWeight=process.effSFLepton2DEventWeight.clone()
    process.effSFElectronEventWeight.particles=cms.InputTag("signalElectrons")
    #process.effSFElectronEventWeight.jets=cms.InputTag("signalTightJets")
    process.effSFElectronEventWeight.sysVar   = cms.string("")
    #process.effSFElectronEventWeight.filename= "TopAnalysis/Configuration/data/EleEffSF2D2012.root"
    
    #process.effSFElectronEventWeightLepID     = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root")
    process.effSFElectronEventWeightLepID     = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/egammaEffi.txt_SF2D_ID_Loose.root")
    
    process.effSFElectronEventWeightLepID.histoname="EGamma_SF2D"  #pt>120 is missing, solved via WeightEventAnalyzer
    
    process.effSFElectronEventWeightUpLepID  =process.effSFElectronEventWeightLepID.clone(sysVar = "combinedEffSFNormUpStat")
    process.effSFElectronEventWeightDownLepID=process.effSFElectronEventWeightLepID.clone(sysVar = "combinedEffSFNormDownStat")
    
	##process.effSFElectronEventWeightIso     = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/eleIsolation_SF.root")
	##process.effSFElectronEventWeightIso.histoname="IsolationSF"  #pt>120 is missing, solved via WeightEventAnalyzer
	#process.effSFElectronEventWeightIso     = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/Electron_IdIso0p15_eff.root")
	#process.effSFElectronEventWeightIso.histoname="IsolationSF"  #pt>120 is missing, solved via WeightEventAnalyzer
	
	#process.effSFElectronEventWeightUpIso  =process.effSFElectronEventWeightIso.clone(sysVar = "combinedEffSFNormUpStat")
	#process.effSFElectronEventWeightDownIso=process.effSFElectronEventWeightIso.clone(sysVar = "combinedEffSFNormDownStat")
	
	##process.effSFElectronEventWeightTrigger    = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/eleTrig_SF.root")
	##process.effSFElectronEventWeightTrigger.histoname="h_eleTrig_SF"  #pt>120 is missing, solved via WeightEventAnalyzer

	##process.effSFElectronEventWeightTrigger    = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/eleRECO.txt_SF2D.root")
	#process.effSFElectronEventWeightTrigger    = process.effSFElectronEventWeight.clone(filename="TopMass/data/LeptonSFs/Electron_Ele25eta2p1WPTight_eff.root")
	
	#process.effSFElectronEventWeightTrigger.histoname="EGamma_SF2D"  #pt>120 is missing, solved via WeightEventAnalyzer
	#process.effSFElectronEventWeightTrigger.additionalSystErr=0  #+ 0.5% for single muon triggers to be applyied on top of statistical errors from root/pkl files
	#process.effSFElectronEventWeightTrigger.verbose=cms.int32(0) 
	
	#process.effSFElectronEventWeightUpTrigger  =process.effSFElectronEventWeightTrigger.clone(sysVar = "combinedEffSFNormUpStat")
	#process.effSFElectronEventWeightDownTrigger=process.effSFElectronEventWeightTrigger.clone(sysVar = "combinedEffSFNormDownStat")
	
	#process.effSFElectronEventWeight.filename= "TopMass/data/LeptonSFs/eleTrig_SF.root" #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHbbRun2ReferenceAnalysis_76XTransition
	#process.effSFElectronEventWeight.histoname= "h_eleTrig_SF"  #pt>120 is missing, solved via WeightEventAnalyzer
	#process.effSFElectronEventWeight.additionalSystErr= 0  #pt>120 is missing, solved via WeightEventAnalyzer
	#process.effSFElectronEventWeight.verbose=cms.int32(0) 
	
	#process.effSFElectronEventWeightUp  =process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormUpStat")
	#process.effSFElectronEventWeightDown=process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormDownStat")

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = outputNTupel
)

# register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = ""

#process.MessageLogger = cms.Service("MessageLogger")



process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
   maxEventsToPrint = cms.untracked.int32(20),
   printVertex = cms.untracked.bool(False),
   src = cms.InputTag("genParticles")
)

