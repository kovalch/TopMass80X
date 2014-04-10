import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os


####################################################################
## Global job options

REPORTEVERY = 100
WANTSUMMARY = True


####################################################################
## Define the process

process = cms.Process("topDileptonNtuple")
#SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )



####################################################################
## Set up command line options

options = VarParsing.VarParsing ('standard')
options.register('runOnMC', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "decide to run on MC or data")
options.register('runOnAOD', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "run on AOD")
options.register('globalTag', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which globalTag should be used")
options.register('mode', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which type of analysis to run")
options.register('samplename', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which sample to run over")
options.register('inputScript', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "python file with input source")
options.register('outputFile', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "root output file")
options.register('systematicsName', 'Nominal', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "type of systematics")
options.register('json', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "limit to certain lumis")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('includePDFWeights', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "include the PDF weights *slow!!!*")

## Get and parse the command line arguments
if( hasattr(sys, "argv") ):
    for args in sys.argv :
        arg = args.split(',')
        for val in arg:
            val = val.split('=')
            if(len(val)==2):
                setattr(options,val[0], val[1])



####################################################################
## Set up samplename

if options.samplename == '':
    print 'cannot run without specifying a samplename'
    exit(8888)

if options.samplename == 'data':
    options.runOnMC = False



####################################################################
## Define input

if options.inputScript != '':
    process.load(options.inputScript)
else:
    print 'need an input script'
    exit(8889)

print "max events: ", options.maxEvents
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

if options.skipEvents > 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)

# Limit to json file (if passed as parameter)
if options.json != '':
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = options.json).getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)



####################################################################
## Create output path

if options.outputFile == '':
    fn = options.mode + '_test.root'
else:
    fn = options.outputFile
print 'Using output file ' + fn

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(fn)
)



####################################################################
## Configure message logger

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = REPORTEVERY

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(WANTSUMMARY)
)



####################################################################
## Geometry and Detector Conditions

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if options.globalTag != '':
    print "Setting global tag to the command-line value"
    process.GlobalTag.globaltag = cms.string( options.globalTag )
else:
    print "Determine global tag automatically"
    if options.runOnMC:
        process.GlobalTag.globaltag = cms.string('START53_V27::All')
    else:
        process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')

print "Using global tag: ", process.GlobalTag.globaltag

process.load("Configuration.StandardSequences.MagneticField_cff")



####################################################################
## HCAL Noise filter

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)



####################################################################
## Beam scraping filter

process.scrapingFilter = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn     = cms.untracked.bool(False),
    numtrack    = cms.untracked.uint32(10),
    thresh      = cms.untracked.double(0.25)
    )



####################################################################
## ECAL laser correction filter

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")



####################################################################
## Primary vertex filtering

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )



####################################################################
## Trigger filtering

# Get the central diLepton trigger lists, and set up filter
from TopAnalysis.TopFilter.sequences.diLeptonTriggers_cff import *
process.load("TopAnalysis.TopFilter.filters.TriggerFilter_cfi")
process.filterTrigger.TriggerResults = cms.InputTag('TriggerResults','','HLT')
process.filterTrigger.printTriggers = False
if options.mode == 'mumu':
    process.filterTrigger.hltPaths  = mumuTriggers
elif options.mode == 'emu':
    process.filterTrigger.hltPaths  = emuTriggers
elif options.mode == 'ee':
    process.filterTrigger.hltPaths  = eeTriggers
else:
    process.filterTrigger.hltPaths = eeTriggers + emuTriggers + mumuTriggers
    
#print "Printing triggers: ", process.filterTrigger.printTriggers



####################################################################
## Jet corrections

if options.runOnMC:
    jetCorr = ('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
else:
    jetCorr = ('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute', 'L2L3Residual'])



####################################################################
## PF2PAT sequence

process.load("Configuration.EventContent.EventContent_cff")
process.out = cms.OutputModule("PoolOutputModule",
    process.FEVTEventContent,
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('RECO')),
    fileName = cms.untracked.string("eh.root"),
)

process.load("PhysicsTools.PatAlgos.patSequences_cff")


pfpostfix = "PFlow"

from PhysicsTools.PatAlgos.tools.pfTools import *

usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=options.runOnMC, postfix=pfpostfix, jetCorrections=jetCorr, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True) 

from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *

# Produce pat trigger content
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")



####################################################################
## Set up selections for PF2PAT & PAT objects: Electrons

# MVA ID
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.electronMvaIdSequence = cms.Sequence(process.mvaTrigV0 + process.mvaNonTrigV0)
getattr(process,'patElectrons'+pfpostfix).electronIDSources.mvaTrigV0 = cms.InputTag("mvaTrigV0")
getattr(process,'patElectrons'+pfpostfix).electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
getattr(process, 'patPF2PATSequence'+pfpostfix).replace(getattr(process,'patElectrons'+pfpostfix),
                                                process.electronMvaIdSequence *
                                                getattr(process,'patElectrons'+pfpostfix)
                                                )


#default is:
#process.pfSelectedElectrons.cut = 'pt > 5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'
#here a low pt cut needs to be used!
getattr(process, 'pfSelectedElectrons'+pfpostfix).cut = 'pt > 5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfHits <= 0'


# Switch isolation cone to 0.3 and set cut to 0.15
getattr(process, 'pfIsolatedElectrons'+pfpostfix).doDeltaBetaCorrection = True   # not really a 'deltaBeta' correction, but it serves
getattr(process, 'pfIsolatedElectrons'+pfpostfix).deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId"+pfpostfix)
getattr(process, 'pfIsolatedElectrons'+pfpostfix).isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"+pfpostfix))
getattr(process, 'pfIsolatedElectrons'+pfpostfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"+pfpostfix), cms.InputTag("elPFIsoValueGamma03PFId"+pfpostfix))
getattr(process, 'pfIsolatedElectrons'+pfpostfix).isolationCut = 0.15


getattr(process, 'patElectrons'+pfpostfix).isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId"+pfpostfix),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId"+pfpostfix),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId"+pfpostfix),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId"+pfpostfix),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId"+pfpostfix)
)


getattr(process, 'selectedPatElectrons'+pfpostfix).cut = 'electronID("mvaTrigV0") > 0.5 && passConversionVeto'
            #cant do this on python level :-(
            #' && abs(gsfTrack().dxy(vertex_.position())) < 0.04'


process.selectedPatElectronsAfterScaling = selectedPatElectrons.clone(
    src = 'scaledJetEnergy:selectedPatElectrons'+pfpostfix,
    cut = 'pt > 20 && abs(eta) < 2.5'
)



####################################################################
## Set up selections for PF2PAT & PAT objects: Muons

getattr(process, 'pfSelectedMuons'+pfpostfix).cut = 'pt > 5' \
                              '&& muonRef.isNonnull()' \
                              '&& (muonRef.isGlobalMuon() || muonRef.isTrackerMuon())'

# Switch isolation cone to 0.3 and set cut to 0.15
getattr(process, 'pfIsolatedMuons'+pfpostfix).doDeltaBetaCorrection = True
getattr(process, 'pfIsolatedMuons'+pfpostfix).deltaBetaIsolationValueMap = cms.InputTag("muPFIsoValuePU03"+pfpostfix, "", "")
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationValueMapsCharged = [cms.InputTag("muPFIsoValueCharged03"+pfpostfix)]
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationValueMapsNeutral = [cms.InputTag("muPFIsoValueNeutral03"+pfpostfix), cms.InputTag("muPFIsoValueGamma03"+pfpostfix)]
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationCut = 0.15

getattr(process,'patMuons'+pfpostfix).isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral03"+pfpostfix),
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll03"+pfpostfix),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU03"+pfpostfix),
        pfPhotons = cms.InputTag("muPFIsoValueGamma03"+pfpostfix),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged03"+pfpostfix)
)


getattr(process,'selectedPatMuons'+pfpostfix).cut = 'isPFMuon && pt > 20 && abs(eta) < 2.4'



####################################################################
## Set up selections for PF2PAT & PAT objects: Jets

getattr(process,'selectedPatJets'+pfpostfix).cut = 'abs(eta)<5.4'


####################################################################
##  MVA met

process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
process.calibratedAK5PFJetsForPFMEtMVA.src = 'pfJets' +pfpostfix
if options.runOnMC:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors= cms.vstring("ak5PFL1FastL2L3")
else:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors= cms.vstring("ak5PFL1FastL2L3Residual")
process.pfMEtMVA.srcUncorrJets='pfJets' +pfpostfix
process.pfMEtMVA.srcVertices = 'goodOfflinePrimaryVertices'
process.pfMEtMVA.inputFileNames = cms.PSet(
    U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_53_June2013_type1.root'),
    DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
    CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
    CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root')
    )
process.pfMEtMVA.srcLeptons=cms.VInputTag("isomuons","isoelectrons","isotaus") #should be adapted to analysis selection..
process.pfMEtMVA.srcRho = cms.InputTag("kt6PFJets","rho","RECO")
process.patMEtMVA= getattr(process,'patMETs'+pfpostfix).clone()
process.patMEtMVA.metSource = 'pfMEtMVA'


####################################################################
## Phi correction of the PFMET
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data = cms.PSet( # CV: ReReco data + Summer'13 JEC
    px = cms.string("+4.83642e-02 + 2.48870e-01*Nvtx"),
    py = cms.string("-1.50135e-01 - 8.27917e-02*Nvtx")
)

process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc = cms.PSet( # CV: Summer'12 MC + Summer'13 JEC
    px = cms.string("+1.62861e-01 - 2.38517e-02*Nvtx"),
    py = cms.string("+3.60860e-01 - 1.30335e-01*Nvtx")
)

process.pfMEtSysShiftCorr.src = cms.InputTag('pfMet')

if options.runOnMC:
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
else:
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data



####################################################################
## Basic debugging analyzer

#process.load("TopAnalysis.TopAnalyzer.CheckDiLeptonAnalyzer_cfi")
#process.analyzeDiLepton.electrons = 'fullySelectedPatElectronsCiC'
#process.analyzeDiLepton.muons = 'fullySelectedPatMuons'



####################################################################
## Set up sample-specific flags for individual treatment in nTuple production

zGenInfo = False
zproducer = False
topfilter = False
signal = False
higgsSignal = False
alsoViaTau = False
ttbarV = False

if options.samplename == 'ttbarsignal':
    topfilter = True
    signal = True
    viaTau = False
elif options.samplename == 'ttbarsignalviatau':
    topfilter = True
    signal = True
    viaTau = True
elif options.samplename == 'ttbarsignalplustau':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
elif options.samplename == 'ttbarbg':
    topfilter = True
elif options.samplename == 'dy1050' or options.samplename == 'dy50inf':
    zproducer = True
    zGenInfo = True
elif options.samplename == 'ttbarhiggstobbbar' or options.samplename == 'ttbarhiggsinclusive':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
    higgsSignal = True
elif options.samplename == 'gghiggstozzto4l' or options.samplename == 'vbfhiggstozzto4l':
    zGenInfo = True
    higgsSignal = True
elif options.samplename == 'ttbarw' or options.samplename == 'ttbarz':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
    ttbarV = True
elif options.samplename in ['data', 'singletop', 'singleantitop','ww',
        'wz','zz','wjets',
        'qcdmu15','qcdem2030','qcdem3080','qcdem80170',
        'qcdbcem2030','qcdbcem3080','qcdbcem80170',
        'zzz','wwz','www','ttww','ttg','wwg']:
    #no special treatment needed, put here to avoid typos
    pass
else:
    print "Error: Unknown samplename!"
    exit(8)


if signal:
    process.goodOfflinePrimaryVertices.filter = cms.bool(False)

####################################################################
## Define which collections (including which corrections) to be used in nTuple

isolatedMuonCollection = "selectedPatMuons"+pfpostfix

#isolatedElecCollection = "selectedPatElectrons"
isolatedElecCollection = "selectedPatElectronsAfterScaling"

jetCollection = "hardJets"

jetForMETCollection = "scaledJetEnergy:selectedPatJets"+pfpostfix

metCollection = "scaledJetEnergy:patMETs"+pfpostfix

genJetCollection = "ak5GenJetsPlusHadron"

genLevelBJetProducerInput = "produceGenLevelBJets"

genHFBHadronMatcherInput = "matchGenHFBHadronJets"
genHFCHadronMatcherInput = "matchGenHFCHadronJets"

# Lepton collection used for kinematic reconstruction (has further selections, and can thus deviate from KinReco in analysis)
finalLeptons = 'filterDiLeptonMassQCDveto'



####################################################################
## Separation of ttbar samples in dileptonic and other decays

if topfilter:
    process.load("TopAnalysis.TopFilter.filters.GeneratorTopFilter_cfi")
    process.generatorTopFilter.rejectNonBottomDecaysOfTops = False
    if higgsSignal or ttbarV:
        process.generatorTopFilter.invert_selection = True
        process.generatorTopFilter.channels = ["none"] #empty array would use some defaults
    else:
        all = ['ElectronElectron', 'ElectronElectronViaTau', 
               'MuonMuon', 'MuonMuonViaTau', 
               'ElectronMuon', 'ElectronMuonViaTau']
        if signal:
                process.generatorTopFilter.invert_selection = False
                if viaTau:
                        process.generatorTopFilter.channels = ['ElectronElectronViaTau', 'MuonMuonViaTau', 'ElectronMuonViaTau']
                elif alsoViaTau:
                        process.generatorTopFilter.channels = all
                else:
                        process.generatorTopFilter.channels = ['ElectronElectron', 'ElectronMuon', 'MuonMuon']
        else:
                process.generatorTopFilter.channels = all
                process.generatorTopFilter.invert_selection = True



####################################################################
## Build Jet Collections

process.load("TopAnalysis.TopUtils.JetEnergyScale_cfi")

process.load("TopAnalysis.TopFilter.filters.JetIdFunctorFilter_cfi")
process.goodIdJets.jets    = cms.InputTag("scaledJetEnergy:selectedPatJets"+pfpostfix)
process.goodIdJets.jetType = cms.string('PF')
process.goodIdJets.version = cms.string('FIRSTDATA')
process.goodIdJets.quality = cms.string('LOOSE')

process.hardJets = selectedPatJets.clone(src = 'goodIdJets', cut = 'pt > 5 & abs(eta) < 2.4') 

# Additional properties for jets like jet charges
process.load("TopAnalysis.HiggsUtils.producers.JetPropertiesProducer_cfi")
process.jetProperties.src = jetCollection

process.buildJets = cms.Sequence(
            process.scaledJetEnergy *
            process.selectedPatElectronsAfterScaling *
            process.goodIdJets * 
            process.hardJets *
            process.jetProperties
            )



####################################################################
## Filter on events containing dilepton system of opposite charge and above m(ll) > 12 GeV

##### WARNING: This tool selects lepton pair based on highest pt-sum of the leptons. Since opposite charge is required, it will form the combination of pt-leading lepton and pt-leading antilepton.
##### WARNING: So it might have strange side effects in selections of events with >=3 leptons ?!
##### WARNING: It will for sure have strange side effects in case of additional lepton selections on nTuple level...
##### WARNING: We should replace it by a simple dilepton selector checking for any acceptable combination, not making any specific choice
from TopAnalysis.TopFilter.filters.DiLeptonFilter_cfi import *
process.filterOppositeCharge = filterLeptonPair.clone(
    electrons    = isolatedElecCollection,
    muons        = isolatedMuonCollection,
    Cut          = (0.,0.),
    filterCharge = -1,
)


##### WARNING: This tool uses the lepton pair from the DiLeptonFilter_cfi, based on highest pt-sum of the leptons.
from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
process.filterChannel = countPatLeptons.clone()
process.filterChannel.electronSource    = 'filterOppositeCharge'
process.filterChannel.muonSource        = 'filterOppositeCharge'
process.filterChannel.minNumber         = 2
process.filterChannel.countTaus         = False
if options.mode == 'ee':
    process.filterChannel.countElectrons    = True
    process.filterChannel.countMuons        = False
elif options.mode == 'mumu':
    process.filterChannel.countElectrons    = False
    process.filterChannel.countMuons        = True
elif options.mode == 'emu':
    process.filterChannel.minNumber         = 1
    process.filterChannel1 = process.filterChannel.clone()
    process.filterChannel2 = process.filterChannel1.clone()
    process.filterChannel1.countElectrons    = True
    process.filterChannel1.countMuons        = False
    process.filterChannel2.countElectrons    = False
    process.filterChannel2.countMuons        = True
    process.filterChannel = cms.Sequence(process.filterChannel1 * process.filterChannel2)
else:
    process.filterChannel.countElectrons    = True
    process.filterChannel.countMuons        = True


##### WARNING: This tool uses the lepton pair from the DiLeptonFilter_cfi, based on highest pt-sum of the leptons.
process.filterDiLeptonMassQCDveto           = filterLeptonPair.clone()
process.filterDiLeptonMassQCDveto.muons     = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.electrons = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.Cut       = (0.,12.)



####################################################################
## Write Ntuple

from TopAnalysis.TopAnalyzer.NTupleWriter_cfi import writeNTuple
writeNTuple.sampleName = options.samplename
writeNTuple.channelName = options.mode
writeNTuple.systematicsName = options.systematicsName
writeNTuple.isMC = options.runOnMC
writeNTuple.isTtBarSample = signal
writeNTuple.isHiggsSample = higgsSignal
writeNTuple.isZSample = zGenInfo
writeNTuple.includePDFWeights = options.includePDFWeights
writeNTuple.pdfWeights = "pdfWeights:cteq66"
writeNTuple.includeZdecay = zproducer
writeNTuple.saveHadronMothers = False



process.writeNTuple = writeNTuple.clone(
    muons = isolatedMuonCollection,
    elecs = isolatedElecCollection,
    jets = jetCollection,
    met = metCollection,
    mvamet = "patMEtMVA",
    genMET = "genMetTrue",
    genJets = genJetCollection,

    BHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "BHadJetIndex"),
    AntiBHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "AntiBHadJetIndex"),
    BHadrons = cms.InputTag(genLevelBJetProducerInput, "BHadrons"),
    AntiBHadrons = cms.InputTag(genLevelBJetProducerInput, "AntiBHadrons"),
    BHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "BHadronFromTopB"),
    AntiBHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronFromTopB"),
    BHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "BHadronVsJet"),
    AntiBHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronVsJet"),
    genBHadPlusMothers = cms.InputTag(genHFBHadronMatcherInput,"genBHadPlusMothers"),
    genBHadPlusMothersIndices = cms.InputTag(genHFBHadronMatcherInput,"genBHadPlusMothersIndices"),
    genBHadIndex = cms.InputTag(genHFBHadronMatcherInput,"genBHadIndex"),
    genBHadFlavour = cms.InputTag(genHFBHadronMatcherInput,"genBHadFlavour"),
    genBHadJetIndex = cms.InputTag(genHFBHadronMatcherInput,"genBHadJetIndex"),
    genBHadFromTopWeakDecay = cms.InputTag(genHFBHadronMatcherInput,"genBHadFromTopWeakDecay"),
    genBHadLeptonIndex = cms.InputTag(genHFBHadronMatcherInput,"genBHadLeptonIndex"),
    genBHadLeptonHadronIndex = cms.InputTag(genHFBHadronMatcherInput,"genBHadLeptonHadronIndex"),
    genBHadLeptonViaTau = cms.InputTag(genHFBHadronMatcherInput,"genBHadLeptonViaTau"),

    genCHadPlusMothers = cms.InputTag(genHFCHadronMatcherInput,"genCHadPlusMothers"),
    genCHadPlusMothersIndices = cms.InputTag(genHFCHadronMatcherInput,"genCHadPlusMothersIndices"),
    genCHadIndex = cms.InputTag(genHFCHadronMatcherInput,"genCHadIndex"),
    genCHadFlavour = cms.InputTag(genHFCHadronMatcherInput,"genCHadFlavour"),
    genCHadJetIndex = cms.InputTag(genHFCHadronMatcherInput,"genCHadJetIndex"),
    genCHadFromTopWeakDecay = cms.InputTag(genHFCHadronMatcherInput,"genCHadFromTopWeakDecay"),
    genCHadLeptonIndex = cms.InputTag(genHFCHadronMatcherInput,"genCHadLeptonIndex"),
    genCHadLeptonHadronIndex = cms.InputTag(genHFCHadronMatcherInput,"genCHadLeptonHadronIndex"),
    genCHadLeptonViaTau = cms.InputTag(genHFCHadronMatcherInput,"genCHadLeptonViaTau"), 
)
process.writeNTuple.jetsForMET    = cms.InputTag("scaledJetEnergy:selectedPatJets"+pfpostfix)
process.writeNTuple.jetsForMETuncorr    = cms.InputTag("selectedPatJets"+pfpostfix)



####################################################################
## Include PDF weights for systematic signal samples

if options.includePDFWeights:
    if not signal:
        print "PDF variations only supported for the signal"
        exit(5615)
    process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                # Fix POWHEG if buggy (this PDF set will also appear on output,
                # so only two more PDF sets can be added in PdfSetNames if not "")
                #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                GenTag = cms.untracked.InputTag("genParticles"),
                PdfInfoTag = cms.untracked.InputTag("generator"),
                PdfSetNames = cms.untracked.vstring(
                        "cteq66.LHgrid"
                        #, "MRST2006nnlo.LHgrid"
                        #, "NNPDF10_100.LHgrid"
                        #"cteq6mE.LHgrid"
                        # ,"cteq6m.LHpdf"
                        #"cteq6m.LHpdf"
                ))
else:
    process.pdfWeights = cms.Sequence()



####################################################################
## Kinematic reconstruction

# std sequence to produce the ttFullLepEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff")
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff import *
if not signal:
    removeTtFullLepHypGenMatch(process)

setForAllTtFullLepHypotheses(process,"muons"    ,finalLeptons)
setForAllTtFullLepHypotheses(process,"electrons",finalLeptons)
# WARNING! The jet.pt > 30 cut is currently hardcoded in the NTupleWriter.cc file
# adding a collections like
#     process.jetsForKinReco = process.hardJets.clone(src = 'hardJets', cut = 'pt > 30')
# will cause problems because the selection of the "best" solution is hardcoded!!!!!
#setForAllTtFullLepHypotheses(process,"jets"     ,'jetsForKinReco')
setForAllTtFullLepHypotheses(process,"jets"     ,jetCollection)
setForAllTtFullLepHypotheses(process,"mets"     ,metCollection)
if options.runOnMC:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L3Absolute")
    print "L3Absolute"
else:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L2L3Residual")
    print "L2L3Residual"
setForAllTtFullLepHypotheses(process,"maxNJets",-1)

process.kinSolutionTtFullLepEventHypothesis.maxNComb = -1
process.kinSolutionTtFullLepEventHypothesis.searchWrongCharge = True
process.kinSolutionTtFullLepEventHypothesis.tmassbegin = 100.0
process.kinSolutionTtFullLepEventHypothesis.tmassend   = 300.0
process.kinSolutionTtFullLepEventHypothesis.neutrino_parameters = (30.641, 57.941, 22.344, 57.533, 22.232)
#according to our MC 8 TeV values are:
#nu    mpv 40.567 sigma = 16.876
#nubar mpv 40.639 sigma = 17.021

#process.kinSolutionTtFullLepEventHypothesis.mumuChannel = False
#process.kinSolutionTtFullLepEventHypothesis.eeChannel = False
#process.kinSolutionTtFullLepEventHypothesis.emuChannel = True
process.ttFullLepEvent.decayChannel1 = cms.int32(1)
process.ttFullLepEvent.decayChannel2 = cms.int32(2)

process.kinSolutionTtFullLepEventHypothesis.mumuChannel = True
process.kinSolutionTtFullLepEventHypothesis.emuChannel  = True
process.kinSolutionTtFullLepEventHypothesis.eeChannel = True

#process.ttFullLepEvent.decayChannel1 = cms.int32(1)
#process.ttFullLepEvent.decayChannel2 = cms.int32(2)


#if options.mode == 'mumu':
    #process.kinSolutionTtFullLepEventHypothesis.mumuChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(2)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(2)
#elif options.mode == 'emu':
    #process.kinSolutionTtFullLepEventHypothesis.emuChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(2)
#elif options.mode == 'ee':
    #process.kinSolutionTtFullLepEventHypothesis.eeChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(1)



####################################################################
## Sample-specific sequences

if zproducer:
    process.load("TopAnalysis.TopUtils.ZDecayProducer_cfi")
    process.zsequence = cms.Sequence(process.ZDecayProducer)
else:
    process.zsequence = cms.Sequence()

if zGenInfo:
    process.load("TopAnalysis.HiggsUtils.producers.GenZDecay_cfi")
    process.zGenSequence = cms.Sequence(process.genZDecay)
else:
    process.zGenSequence = cms.Sequence()

if topfilter:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    #process.load("TopAnalysis.TopUtils.HadronLevelBJetProducer_cfi")

    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # supplies PDG ID to real name resolution of MC particles, necessary for GenLevelBJetProducer
    process.load("TopAnalysis.TopUtils.GenLevelBJetProducer_cfi")
    process.produceGenLevelBJets.deltaR = 5.0
    process.produceGenLevelBJets.noBBbarResonances = True

    from TopAnalysis.TopUtils.GenHFHadronMatcher_cff import matchGenHFHadronJets
    process.matchGenHFBHadronJets = matchGenHFHadronJets.clone(
        flavour = 5,
        noBBbarResonances = True
    )
    process.matchGenHFCHadronJets = matchGenHFHadronJets.clone(
        flavour = 4,
        noBBbarResonances = True
    )

    process.load("TopAnalysis.TopUtils.sequences.improvedJetHadronQuarkMatching_cff")

    process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
    if signal:
        process.topsequence = cms.Sequence(
            process.makeGenEvt *
            process.improvedJetHadronQuarkMatchingSequence *
            process.generatorTopFilter *
            process.produceGenLevelBJets *
            process.matchGenHFBHadronJets *
            process.matchGenHFCHadronJets)
    else:
        process.topsequence = cms.Sequence(
            process.makeGenEvt *
            process.generatorTopFilter)

else:
    process.topsequence = cms.Sequence()


if higgsSignal:
    process.load("TopAnalysis.HiggsUtils.filters.GeneratorHiggsFilter_cfi")
    process.generatorHiggsFilter.channels = ["none"]
    process.generatorHiggsFilter.invert_selection = True
    process.load("TopAnalysis.HiggsUtils.sequences.higgsGenEvent_cff")
    process.decaySubsetHiggs.fillMode = "kME" # Status3, use kStable for Status2
    process.higgssequence = cms.Sequence(
        process.makeGenEvtHiggs *
        process.generatorHiggsFilter
    )
else:
    process.higgssequence = cms.Sequence()

if signal or higgsSignal or zGenInfo:
    process.ntupleInRecoSeq = cms.Sequence()
else:
    process.ntupleInRecoSeq = cms.Sequence(process.zsequence * process.writeNTuple)
    


####################################################################
## Remove all the tau stuff

from PhysicsTools.PatAlgos.tools.coreTools import removeSpecificPATObjects
removeSpecificPATObjects( process
                        , names = ['Taus', 'Photons']
                        , outputModules = []
                        , postfix = pfpostfix
                        )

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatTaus'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatTaus'+pfpostfix))

## removal of unnecessary modules
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patCandidateSummary'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatCandidateSummary'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatElectrons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatMuons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatLeptons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatJets'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatPFParticles'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfPhotonSequence'+pfpostfix))

massSearchReplaceAnyInputTag(getattr(process,'patPF2PATSequence'+pfpostfix),'pfNoTau'+pfpostfix,'pfJets'+pfpostfix)

#Set to true the access to some specific information in PAT jets when doing the b-tagging. 
applyPostfix(process, "patJets", pfpostfix).addTagInfos = True


####################################################################
## Paths, one with preselection, one without for signal samples

process.p = cms.Path(
    process.goodOfflinePrimaryVertices *
    getattr(process,'patPF2PATSequence'+pfpostfix) *
    process.pfMEtSysShiftCorrSequence     *
    process.buildJets                     *
    process.filterOppositeCharge          *
    process.filterChannel                 *
    process.filterDiLeptonMassQCDveto     *
    process.makeTtFullLepEvent            *
    process.ntupleInRecoSeq
)

if signal or higgsSignal or zGenInfo:
    process.pNtuple = cms.Path(
        process.goodOfflinePrimaryVertices *
        getattr(process,'patPF2PATSequence'+pfpostfix) *
        process.pfMEtSysShiftCorrSequence *
        process.buildJets *
        process.zsequence *
        process.writeNTuple
        )


getattr(process, 'p').replace(process.filterDiLeptonMassQCDveto,
                            process.filterDiLeptonMassQCDveto * process.pfMEtMVAsequence * process.patMEtMVA)

if signal or higgsSignal or zGenInfo:
    getattr(process, 'pNtuple').replace(process.zsequence,
                                        process.zsequence * process.pfMEtMVAsequence * process.patMEtMVA)

####################################################################
## Prepend PF2PAT

from TopAnalysis.TopAnalyzer.CountEventAnalyzer_cfi import countEvents
process.EventsBeforeSelection = countEvents.clone()
process.EventsBeforeSelection.includePDFWeights = options.includePDFWeights
process.EventsBeforeSelection.pdfWeights = "pdfWeights:cteq66"


pathnames = process.paths_().keys()
print 'prepending trigger sequence to paths:', pathnames
for pathname in pathnames:
    getattr(process, pathname).insert(0, cms.Sequence(
        process.pdfWeights *
        process.EventsBeforeSelection *
        process.topsequence *
        process.higgssequence *
        process.zGenSequence *
        process.filterTrigger
        ))
if signal or higgsSignal or zGenInfo:
    process.pNtuple.remove(process.filterTrigger)

process.scaledJetEnergy.inputElectrons       = "selectedPatElectrons"+pfpostfix
process.scaledJetEnergy.inputJets            = "selectedPatJets"+pfpostfix
process.scaledJetEnergy.inputMETs            = "patMETs"+pfpostfix
process.scaledJetEnergy.JECUncSrcFile        = cms.FileInPath("TopAnalysis/Configuration/analysis/common/data/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt")
process.scaledJetEnergy.scaleType = "abs"   #abs = 1, jes:up, jes:down

if options.runOnMC:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, 0.5, 0.5, 1.1, 1.1, 1.7, 1.7, 2.3, 2.3, 5.4)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.052, 1.057, 1.096, 1.134, 1.288) # JER standard

    #please change this on the top where the defaults for the VarParsing are given
    if options.systematicsName == "JES_UP":
        process.scaledJetEnergy.scaleType = "jes:up"
    if options.systematicsName == "JES_DOWN":
        process.scaledJetEnergy.scaleType = "jes:down"
    if options.systematicsName == "JER_UP":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(1.115, 1.114, 1.161, 1.228, 1.488)
    if options.systematicsName == "JER_DOWN":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(0.990, 1.001, 1.032, 1.042, 1.089)
else:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, -1)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.0) # JER standard
    for pathname in pathnames:
        getattr(process, pathname).replace(process.goodOfflinePrimaryVertices,
                                           process.HBHENoiseFilter * 
                                           process.scrapingFilter * 
                                           process.ecalLaserCorrFilter * 
                                           process.goodOfflinePrimaryVertices)
        


####################################################################
## Particle tree drawer

# see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCandidateModules#ParticleTreeDrawer_Utility
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#             #                      printP4 = cms.untracked.bool(False),
#             #                      printPtEtaPhi = cms.untracked.bool(False),
#             #                      printVertex = cms.untracked.bool(False),
#             #                      printStatus = cms.untracked.bool(False),
#             #                      printIndex = cms.untracked.bool(False),
#             #                      status = cms.untracked.vint32( 3 )
#                                   )
#process.p = cms.Path(process.printTree)
#process.pNtuple = cms.Path()



####################################################################
## Signal catcher for more information on errors

process.load("TopAnalysis.TopUtils.SignalCatcher_cfi")



