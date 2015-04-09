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
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(WANTSUMMARY)
)
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

# Get and parse the command line arguments
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
## Set up sample-specific flags for individual treatment in nTuple production

zGenInfo = False
topfilter = False
topSignal = False
higgsSignal = False
alsoViaTau = False
ttbarZ = False
containsWs = False

if options.samplename == 'ttbarsignal':
    topfilter = True
    topSignal = True
    viaTau = False
    containsWs = True
elif options.samplename == 'ttbarsignalviatau':
    topfilter = True
    topSignal = True
    viaTau = True
    containsWs = True
elif options.samplename == 'ttbarsignalplustau':
    topfilter = True
    topSignal = True
    viaTau = False
    alsoViaTau = True
    containsWs = True
elif options.samplename == 'ttbarbg':
    topfilter = True
    containsWs = True
elif options.samplename == 'dy1050' or options.samplename == 'dy50inf':
    zGenInfo = True
elif options.samplename == 'ttbarhiggstobbbar':
    topfilter = True
    topSignal = True
    viaTau = False
    alsoViaTau = True
    higgsSignal = True
elif options.samplename == 'ttbarhiggsinclusive':
    topfilter = True
    topSignal = True
    viaTau = False
    alsoViaTau = True
    higgsSignal = True
    containsWs = True
elif options.samplename == 'gghiggstozzto4l' or options.samplename == 'vbfhiggstozzto4l':
    zGenInfo = True
    higgsSignal = True
elif options.samplename == 'ttbarz':
    zGenInfo = True
    topfilter = True
    topSignal = True
    viaTau = False
    alsoViaTau = True
    ttbarZ = True
    containsWs = True
elif options.samplename in [
        'singletop', 'singleantitop',
        'wjets', 'ww', 'wz',
        'wwz', 'www', 'wwg',
        'ttww', 'ttg', 'ttbarw']:
    containsWs = True
elif options.samplename in [
        'data', 'zz', 'zzz',
        'qcdmu15', 'qcdem2030', 'qcdem3080', 'qcdem80170',
        'qcdbcem2030', 'qcdbcem3080', 'qcdbcem80170']:
    # No special treatment needed, put here to avoid typos
    pass
else:
    print "Error: Unknown samplename!"
    exit(8)

signal = topSignal or higgsSignal or zGenInfo



####################################################################
## Configure message logger

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = REPORTEVERY



####################################################################
## Geometry and Detector Conditions

process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

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



####################################################################
## Configure TFileService

if options.outputFile == '':
    fn = options.mode + '_test.root'
else:
    fn = options.outputFile
print 'Using output file ' + fn

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(fn)
)



####################################################################
## Trigger filtering

if signal:
    process.triggerSequence = cms.Sequence()
else:
    # Get the central diLepton trigger lists, and set up filter
    from TopAnalysis.TopFilter.sequences.diLeptonTriggers_cff import *
    process.load("TopAnalysis.TopFilter.filters.TriggerFilter_cfi")
    process.filterTrigger.TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
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
    process.triggerSequence = cms.Sequence(process.filterTrigger)



####################################################################
## Prefilter sequence

if options.runOnMC:
    process.prefilterSequence = cms.Sequence()
else:
    ## HCAL Noise filter
    process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
    process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
    process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
    process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)
    
    ## Beam scraping filter
    process.scrapingFilter = cms.EDFilter("FilterOutScraping",
        applyfilter = cms.untracked.bool(True),
        debugOn     = cms.untracked.bool(False),
        numtrack    = cms.untracked.uint32(10),
        thresh      = cms.untracked.double(0.25)
    )
    
    ## ECAL laser correction filter
    process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
    
    process.prefilterSequence = cms.Sequence(
        process.HBHENoiseFilter *
        process.scrapingFilter *
        process.ecalLaserCorrFilter
    )



####################################################################
## Primary vertex filtering

selectedPrimaryVertices = 'goodOfflinePrimaryVertices'

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone(minNdof = cms.double(4.0), maxZ = cms.double(24.0)),
    src = cms.InputTag('offlinePrimaryVertices')
    )
if signal:
    process.goodOfflinePrimaryVertices.filter = cms.bool(False)



####################################################################
## Jet energy corrections

if options.runOnMC:
    jetCorrections = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
else:
    jetCorrections = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])



####################################################################
## Full configuration for PF2PAT

pfpostfix = "PFlow"

## Output module for edm files (needed for PAT sequence, even if not used in EndPath)
process.load("Configuration.EventContent.EventContent_cff")
process.out = cms.OutputModule("PoolOutputModule",
    process.FEVTEventContent,
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('RECO')),
    fileName = cms.untracked.string("eh.root"),
)

## PF2PAT sequence
# Parameter checkClosestZVertex = False needs to be set to False when using PF Jets with Charged Hadron Subtraction, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorPFnoPU2012
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=options.runOnMC, postfix=pfpostfix, jetCorrections=jetCorrections, pvCollection=cms.InputTag(selectedPrimaryVertices), typeIMetCorrections=True)
getattr(process, 'pfPileUp'+pfpostfix).checkClosestZVertex = False
process.userPatSequence = getattr(process, 'patPF2PATSequence'+pfpostfix)

## Selections for PF2PAT objects: Electrons
# MVA ID: all IDs except of MVA ID are removed, and this needs to be done in case of electron energy correction application (else code crashes)
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.electronMvaIdSequence = cms.Sequence(process.mvaTrigV0 + process.mvaNonTrigV0)
getattr(process, 'patElectrons'+pfpostfix).electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
)
process.userPatSequence.replace(getattr(process, 'patElectrons'+pfpostfix),
                                process.electronMvaIdSequence * getattr(process, 'patElectrons'+pfpostfix))
# Basic selection for PF2PAT
getattr(process, 'pfSelectedElectrons'+pfpostfix).cut = 'pt>5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfHits<=0'
# Switch isolation cone to 0.3 and set cut to 0.15
getattr(process, 'pfIsolatedElectrons'+pfpostfix).doDeltaBetaCorrection = True   # Not really a 'deltaBeta' correction, but it serves
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

## Selections for PF2PAT objects: Muons
# Basic selection for PF2PAT
getattr(process, 'pfSelectedMuons'+pfpostfix).cut = 'pt>5 && muonRef.isNonnull() && (muonRef.isGlobalMuon() || muonRef.isTrackerMuon())'
# Switch isolation cone to 0.3 and set cut to 0.15
getattr(process, 'pfIsolatedMuons'+pfpostfix).doDeltaBetaCorrection = True
getattr(process, 'pfIsolatedMuons'+pfpostfix).deltaBetaIsolationValueMap = cms.InputTag("muPFIsoValuePU03"+pfpostfix, "", "")
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationValueMapsCharged = [cms.InputTag("muPFIsoValueCharged03"+pfpostfix)]
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationValueMapsNeutral = [cms.InputTag("muPFIsoValueNeutral03"+pfpostfix), cms.InputTag("muPFIsoValueGamma03"+pfpostfix)]
getattr(process, 'pfIsolatedMuons'+pfpostfix).isolationCut = 0.15
getattr(process, 'patMuons'+pfpostfix).isolationValues = cms.PSet(
    pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral03"+pfpostfix),
    pfChargedAll = cms.InputTag("muPFIsoValueChargedAll03"+pfpostfix),
    pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU03"+pfpostfix),
    pfPhotons = cms.InputTag("muPFIsoValueGamma03"+pfpostfix),
    pfChargedHadrons = cms.InputTag("muPFIsoValueCharged03"+pfpostfix)
)

## Selections for PF2PAT objects: Jets
getattr(process, 'patJetCorrFactors'+pfpostfix).rho = cms.InputTag("kt6PFJets", "rho", "RECO")
# Set addTagInfos to true to access b-tagging information in PAT jets
from PhysicsTools.PatAlgos.tools.helpers import applyPostfix
applyPostfix(process, "patJets", pfpostfix).addTagInfos = True



####################################################################
## Electron energy scale and resolution corrections

process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")
process.eleRegressionEnergy.inputElectronsTag = 'gsfElectrons'
process.eleRegressionEnergy.inputCollectionType = 0
process.eleRegressionEnergy.useRecHitCollections = True
process.eleRegressionEnergy.produceValueMaps = True

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load("EgammaAnalysis.ElectronTools.calibratedElectrons_cfi")
process.calibratedElectrons.inputElectronsTag = 'gsfElectrons'
process.calibratedElectrons.isMC = options.runOnMC
if options.runOnMC:
    process.calibratedElectrons.inputDataset = "Summer12_LegacyPaper"
process.calibratedElectrons.updateEnergyError = True
process.calibratedElectrons.correctionsType = 2
process.calibratedElectrons.combinationType = 3

process.electronCorrectionSequence = cms.Sequence(
    process.eleRegressionEnergy *
    process.calibratedElectrons
)



####################################################################
## Muon energy scale and resolution corrections
## For safety also done for reco muons (i.e. not pf), in case they are used somewhere as input)

process.load("TopAnalysis.ZTopUtils.correctmuonenergy_cff")
process.correctMuonEnergy.muonSrc = 'pfMuonsFromVertex'+pfpostfix
process.correctMuonEnergy.muonType = "pfMuons"
process.correctMuonEnergy.isMC = options.runOnMC
process.correctMuonEnergy.debug = False
process.userPatSequence.replace(getattr(process,'pfMuonsFromVertex'+pfpostfix),
                                getattr(process,'pfMuonsFromVertex'+pfpostfix) * process.correctMuonEnergy)

process.correctRecoMuonEnergy = process.correctMuonEnergy.clone()
process.correctRecoMuonEnergy.muonSrc = 'muons'
process.correctRecoMuonEnergy.muonType = "recoMuons"



####################################################################
## Phi correction of the PFMET

process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data = cms.PSet(
    px = cms.string("+7.99e-02 + 2.393e-01*Nvtx"),
    py = cms.string("-1.370e-01 - 8.02e-02*Nvtx")
)
process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc = cms.PSet(
    px = cms.string("+0.982e-01 + 1.38e-02*Nvtx"),
    py = cms.string("+1.802e-01 - 1.347e-01*Nvtx")
)

process.pfMEtSysShiftCorr.src = 'pfMET'+pfpostfix
if options.runOnMC:
    process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3"
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
else:
    process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3Residual"
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data

# Type1 correction + phi correction
process.pfMetT1Txy = process.pfType1CorrectedMet.clone()
process.pfMetT1Txy.applyType0Corrections = cms.bool(False)
process.pfMetT1Txy.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1'),
    cms.InputTag('pfMEtSysShiftCorr')
)

process.patpfMetT1Txy = getattr(process, 'patPFMet'+pfpostfix).clone()
process.patpfMetT1Txy.metSource = 'pfMetT1Txy'

# Type0 and type1 correction + phi correction
process.pfMetT0T1Txy = process.pfType1CorrectedMet.clone()
process.pfMetT0T1Txy.applyType0Corrections = cms.bool(True)
process.pfMetT0T1Txy.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1'),
    cms.InputTag('pfMEtSysShiftCorr')
)

process.patpfMetT0T1Txy = getattr(process, 'patPFMet'+pfpostfix).clone()
process.patpfMetT0T1Txy.metSource = 'pfMetT0T1Txy'

# Sequence for full inclusion of phi corrections
process.pfMetPhiCorrectionSequence = cms.Sequence(
    process.pfMEtSysShiftCorrSequence *
    process.producePFMETCorrections *
    #process.pfMetT0T1Txy *
    #process.patpfMetT0T1Txy *
    process.pfMetT1Txy *
    process.patpfMetT1Txy
)

process.userPatSequence.replace(getattr(process, 'patMETs'+pfpostfix),
                                process.pfMetPhiCorrectionSequence * getattr(process, 'patMETs'+pfpostfix))



####################################################################
## Calculate MVA MET

process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
process.calibratedAK5PFJetsForPFMEtMVA.src = 'pfJets'+pfpostfix
if options.runOnMC:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
else:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")

process.pfMEtMVA.srcUncorrJets = 'pfJets'+pfpostfix
process.pfMEtMVA.srcVertices = 'goodOfflinePrimaryVertices'
process.pfMEtMVA.inputFileNames = cms.PSet(
    U = cms.FileInPath('TopAnalysis/Configuration/analysis/common/data/gbrmet_53_Sep2013_type1.root'),
    DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
    CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
    CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root')
    )
process.pfMEtMVA.useType1 = True
process.pfMEtMVA.srcLeptons = cms.VInputTag("isomuons", "isoelectrons", "isotaus") # Should be adapted to analysis selection..
process.pfMEtMVA.srcRho = cms.InputTag("kt6PFJets", "rho", "RECO")

process.patMEtMVA = getattr(process, 'patMETs'+pfpostfix).clone()
process.patMEtMVA.metSource = 'pfMEtMVA'

process.mvaMetSequence = cms.Sequence(process.pfMEtMVAsequence * process.patMEtMVA)



####################################################################
## Phi correction of the MVA MET

process.mvaMEtSysShiftCorrParameters_2012runABCDvsNvtx_data = cms.PSet(
    px = cms.string("+1.376e-01 + 0.503e-01*Nvtx"),
    py = cms.string("-0.936e-01 - 0.695e-01*Nvtx")
)
process.mvaMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc = cms.PSet(
    px = cms.string("-1.612e-01 - 0.124e-01*Nvtx"),
    py = cms.string("+0.139e-01 - 0.303e-01*Nvtx")
)

process.mvaMEtSysShiftCorr = process.pfMEtSysShiftCorr.clone(src = 'pfMEtMVA')    
if options.runOnMC:
    process.mvaMEtSysShiftCorr.parameter = process.mvaMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
else:
    process.mvaMEtSysShiftCorr.parameter = process.mvaMEtSysShiftCorrParameters_2012runABCDvsNvtx_data

process.mvaMetTxy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMEtMVA'),
    applyType0Corrections = cms.bool(False),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('mvaMEtSysShiftCorr')
    ),
    applyType2Corrections = cms.bool(False)
)

process.patMvaMetTxy = getattr(process, 'patPFMet'+pfpostfix).clone()
process.patMvaMetTxy.metSource = 'mvaMetTxy'

process.mvaMEtSysShiftCorrSequence = cms.Sequence(
    process.selectedVerticesForMEtCorr *
    process.mvaMEtSysShiftCorr
    )

# Sequence for full inclusion of phi corrections
process.mvaMetPhiCorrectionSequence = cms.Sequence(
    process.mvaMEtSysShiftCorrSequence *
    process.mvaMetTxy *
    process.patMvaMetTxy
    )

process.mvaMetSequence.replace(process.patMEtMVA,
                               process.mvaMetPhiCorrectionSequence * process.patMEtMVA)



####################################################################
## Remove all unneeded PAT modules (taus and photons are fully removed)

from PhysicsTools.PatAlgos.tools.coreTools import removeSpecificPATObjects
removeSpecificPATObjects(process, names = ['Taus', 'Photons'], outputModules = [], postfix = pfpostfix)

process.userPatSequence.remove(getattr(process, 'selectedPatTaus'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatTaus'+pfpostfix))

process.userPatSequence.remove(getattr(process, 'patPFParticles'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'patCandidateSummary'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'selectedPatPFParticles'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'selectedPatCandidateSummary'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatElectrons'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatMuons'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatLeptons'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatJets'+pfpostfix))
process.userPatSequence.remove(getattr(process, 'countPatPFParticles'+pfpostfix))

process.userPatSequence.remove(getattr(process, 'pfPhotonSequence'+pfpostfix))

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
massSearchReplaceAnyInputTag(process.userPatSequence, 'pfNoTau'+pfpostfix, 'pfJets'+pfpostfix)



####################################################################
## Object preselection for any correction based on final PAT objects

preselectedElectronCollection = 'selectedPatElectrons'+pfpostfix
preselectedMuonCollection = 'selectedPatMuons'+pfpostfix
preselectedJetCollection = 'selectedPatJets'+pfpostfix
preselectedMetCollection = 'patpfMetT1Txy'
preselectedMvaMetCollection = 'patMvaMetTxy'

getattr(process, preselectedElectronCollection).cut = 'electronID("mvaTrigV0")>0.5 && passConversionVeto'
getattr(process, preselectedMuonCollection).cut = 'isPFMuon && pt>20 && abs(eta)<2.4'
getattr(process, preselectedJetCollection).cut = 'abs(eta) < 5.4'



####################################################################
## Configure jet energy resolution, and use corrected electrons, and propagate to MET
# Only JER is applied here, JES is already corrected in PAT sequence
# However, an existing file needs to be provided for "JECUncSrcFile", else it crashes (although results do not depend on it)

process.load("TopAnalysis.TopUtils.JetEnergyScale_cfi")
process.scaledJetEnergy.inputElectrons = preselectedElectronCollection
process.scaledJetEnergy.inputJets = preselectedJetCollection
process.scaledJetEnergy.inputMETs = preselectedMetCollection
process.scaledJetEnergy.JECUncSrcFile = cms.FileInPath("TopAnalysis/Configuration/analysis/common/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt")
process.scaledJetEnergy.scaleType = "abs"   # abs = 1, jes:up, jes:down

if options.runOnMC:
    process.scaledJetEnergy.resolutionEtaRanges = cms.vdouble(0, 0.5, 0.5, 1.1, 1.1, 1.7, 1.7, 2.3, 2.3, 2.8, 2.8, 3.2, 3.2, 5.0)
    process.scaledJetEnergy.resolutionFactors = cms.vdouble(1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056)
else:
    process.scaledJetEnergy.resolutionEtaRanges = cms.vdouble(0, -1)
    process.scaledJetEnergy.resolutionFactors = cms.vdouble(1.0)



####################################################################
## Build Final Collections

# Electrons
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
process.selectedElectrons = selectedPatElectrons.clone(
    src = 'scaledJetEnergy:'+preselectedElectronCollection,
    cut = 'pt>20 && abs(eta)<2.5'
)

# Muons
# -- no additional selection required

# Electrons and muons dxy and dz selection
from TopAnalysis.TopUtils.LeptonVertexSelector_cfi import *
process.leptonVertexSelector = leptonVertexSelector.clone()
process.leptonVertexSelector.electrons = "selectedElectrons"
process.leptonVertexSelector.muons = preselectedMuonCollection
process.leptonVertexSelector.vertices = selectedPrimaryVertices
process.leptonVertexSelector.electronDxyMax = 0.04

# Jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
process.load("TopAnalysis.TopFilter.filters.JetIdFunctorFilter_cfi")
process.goodIdJets.jets = cms.InputTag("scaledJetEnergy:"+preselectedJetCollection)
process.goodIdJets.jetType = cms.string('PF')
process.goodIdJets.version = cms.string('FIRSTDATA')
process.goodIdJets.quality = cms.string('LOOSE')
process.selectedJets = selectedPatJets.clone(
    src = 'goodIdJets',
    cut = 'pt>12 & abs(eta)<2.4'
)

process.finalCollectionsSequence = cms.Sequence(
    process.scaledJetEnergy *
    process.selectedElectrons *
    process.leptonVertexSelector *
    process.goodIdJets * 
    process.selectedJets
)



####################################################################
## Final RECO object collections to be used in nTuple

isolatedElectronCollection = "leptonVertexSelector"

isolatedMuonCollection = "leptonVertexSelector"

jetCollection = "selectedJets"

jetForMetUncorrectedCollection = preselectedJetCollection
jetForMetCollection = "scaledJetEnergy:"+jetForMetUncorrectedCollection

metCollection = "scaledJetEnergy:"+preselectedMetCollection

mvaMetCollection = preselectedMvaMetCollection



####################################################################
## Event preselection based on final input collections

if signal:
    process.preselectionSequence = cms.Sequence()
else:
    # Filter on events containing dilepton system of opposite charge and above m(ll) > 12 GeV
    from TopAnalysis.TopFilter.filters.DileptonPreselection_cfi import *
    process.dileptonPreselection = dileptonPreselection.clone(
        electrons = isolatedElectronCollection,
        muons = isolatedMuonCollection,
        filterCharge = -1,
        filterChannel = options.mode,
        excludeMasses = (-999., 12.),
    )
    process.preselectionSequence = cms.Sequence(process.dileptonPreselection)



####################################################################
## Additional properties based on final input collections

# Properties for jets like jet charges
process.load("TopAnalysis.HiggsUtils.producers.JetPropertiesProducer_cfi")
process.jetProperties.src = jetCollection



####################################################################
## Final GEN object collections to be used in nTuple

genParticleCollection = 'genParticles'
genJetInputParticleCollection = 'genParticles'

genJetCollection = 'ak5GenJetsNoNuNoLepton'

genJetFlavourInfoCollection = 'ak5GenJetFlavourPlusLeptonInfos'

genLevelBJetProducerInput = 'produceGenLevelBJets'

genBHadronMatcherInput = 'matchGenBHadron'
genCHadronMatcherInput = 'matchGenCHadron'



####################################################################
## Form gen jets, and jet flavour info with ghost hadrons and leptons injected
## Details in: PhysicsTools/JetExamples/test/printJetFlavourInfo.cc, PhysicsTools/JetExamples/test/printJetFlavourInfo.py
## and in: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition

if topfilter:
    # Supply PDG ID to real name resolution of MC particles
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

    # Input particles for gen jets (stable gen particles to be used in clustering, excluding electrons, muons and neutrinos from hard interaction)
    from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJets
    process.genParticlesForJetsNoNuNoLepton = genParticlesForJets.clone(
        src = genJetInputParticleCollection,
        excludeResonances = True,
        excludeFromResonancePids = [11, 12, 13, 14, 16],
    )

    # Gen jets
    from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
    process.ak5GenJetsNoNuNoLepton = ak5GenJets.clone(src = "genParticlesForJetsNoNuNoLepton")

    # Ghost particle collection for matching to gen jets (b/c hadrons + leptons)
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(particles = genParticleCollection)

    # Flavour info: jet collection with all associated ghosts
    # For the moment leptons need to be specified explicitely here, until lepton access can be made more generic in miniAOD
    # This is only needed as long as the jetConstituents are not accessible directly in miniAOD, then it should be fixed
    # by using the leptons from the constituents, instead of feeding them as ghosts into the jets
    from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
    process.ak5GenJetFlavourPlusLeptonInfos = ak5JetFlavourInfos.clone(
        jets = genJetCollection,
        leptons = cms.InputTag("selectedHadronsAndPartons", "leptons")
    )

    process.genJetSequence = cms.Sequence(
            process.genParticlesForJetsNoNuNoLepton *
            process.ak5GenJetsNoNuNoLepton *
            process.selectedHadronsAndPartons *
            process.ak5GenJetFlavourPlusLeptonInfos
    )
else:
    process.genJetSequence = cms.Sequence()



####################################################################
## Separation of ttbar samples in dileptonic and other decays

pythia8Sample = False
if options.inputScript.find("_pythia8_") == -1:
    pass
else:
    pythia8Sample = True


if topfilter:
    # FIXME: ttGenEvent is now fixed to work on Pythia8, but the fix is a local backport of CMSSW_7 version
    # FIXME: Need to introduce hack to allow application for Pythia6 samples in Pythia8 mode, to compare with old implementation
    # FIXME: Commented line is old implementation, kept for comparisons (not working with Pythia8)
    #process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    process.load("TopAnalysis.TopUtils.sequences.ttGenEventLocal_cff")
    process.initSubset.src = genParticleCollection
    process.decaySubset.src = genParticleCollection
    if pythia8Sample:
        process.decaySubset.runMode = "Run2"
        process.decaySubset.fillMode = "kStable" # Top before Decay, after Radiation
    else:
        process.decaySubset.runMode = "Run1"
        process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
    
    if not pythia8Sample:
        # FIXME: The switch for Pythia8 is important because of status codes
        process.load("TopAnalysis.TopFilter.filters.GeneratorTopFilter_cfi")
    else:
        process.load("TopAnalysis.TopFilter.filters.GeneratorTopFilter_Pythia8_cfi")
        process.generatorTopFilter.src = genParticleCollection
    process.generatorTopFilter.rejectNonBottomDecaysOfTops = False
    if higgsSignal or ttbarZ:
        process.generatorTopFilter.invert_selection = True
        process.generatorTopFilter.channels = ["none"] # Empty array would use some defaults
    else:
        all = ['ElectronElectron', 'ElectronElectronViaTau', 
               'MuonMuon', 'MuonMuonViaTau', 
               'ElectronMuon', 'ElectronMuonViaTau']
        if topSignal:
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
    
    process.topSplittingSequence = cms.Sequence(
        process.makeGenEvt *
        process.generatorTopFilter
    )
else:
    process.topSplittingSequence = cms.Sequence()



####################################################################
## Sample-specific sequences for generator level information

if zGenInfo:
    process.load("TopAnalysis.HiggsUtils.producers.GenZDecay_cfi")
    process.genZDecay.src = genParticleCollection
    process.zGenSequence = cms.Sequence(process.genZDecay)
else:
    process.zGenSequence = cms.Sequence()

if topSignal:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # Supplies PDG ID to real name resolution of MC particles
    from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenBHadron
    from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenCHadron
    process.matchGenBHadron = matchGenBHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = genJetFlavourInfoCollection,
        onlyJetClusteredHadrons = False,
    )
    process.matchGenCHadron = matchGenCHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = genJetFlavourInfoCollection,
        onlyJetClusteredHadrons = False,
    )
    
    # FIXME: Do we need this, or is this now obsolete? If yes, needs to be adjusted to official GenHFHadronMatcher  
    #process.genParticlesForJetsPlusBHadron.excludeResonances = False
    #process.genParticlesForJetsPlusCHadron.excludeResonances = False
    #process.genParticlesForJetsPlusBCHadron.excludeResonances = False

    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # Supplies PDG ID to real name resolution of MC particles
    process.load("TopAnalysis.TopUtils.GenLevelBJetProducer_cfi")
    process.produceGenLevelBJets.deltaR = 5.0
    process.produceGenLevelBJets.noBBbarResonances = True
    process.produceGenLevelBJets.genJets = genJetCollection

    process.genJetFlavourSequence = cms.Sequence(
        process.matchGenBHadron *
        process.matchGenCHadron *
        process.produceGenLevelBJets
    )
else:
    process.genJetFlavourSequence = cms.Sequence()

if higgsSignal:
    process.load("TopAnalysis.HiggsUtils.sequences.higgsGenEvent_cff")
    process.decaySubsetHiggs.src = genParticleCollection
    process.decaySubsetHiggs.fillMode = "kME" # Status3, use kStable for Status2
    ## FIXME: Workaround for Pythia8, higgsGenEvent not working there
    if not pythia8Sample:
        process.load("TopAnalysis.HiggsUtils.filters.GeneratorHiggsFilter_cfi")
    else:
        process.load("TopAnalysis.HiggsUtils.filters.GeneratorHiggsFilter_Pythia8_cfi")
        process.generatorHiggsFilter.src = genParticleCollection 
    process.generatorHiggsFilter.channels = ["none"]
    process.generatorHiggsFilter.invert_selection = True
    
    process.higgsGenSequence = cms.Sequence(
        process.makeGenEvtHiggs *
        process.generatorHiggsFilter
    )
else:
    process.higgsGenSequence = cms.Sequence()



####################################################################
## W decay modes for MadGraph samples, in order to correct branching ratios

madgraphSample = False
process.madgraphWDecaySequence = cms.Sequence()
if containsWs:
    if options.inputScript.find("_madgraph_" or "_madgraphMLM_" or "_amcatnloFXFX_") == -1:
        pass
    else:
        madgraphSample = True
        process.load("TopAnalysis.TopUtils.MadgraphWDecayProducer_cfi")
        process.madgraphWDecayProducer.src = genParticleCollection
        process.madgraphWDecaySequence = cms.Sequence(process.madgraphWDecayProducer)



####################################################################
## Include PDF weights for systematic signal samples

if options.includePDFWeights:
    if not topSignal:
        print "PDF variations only supported for the signal"
        exit(5615)
    process.pdfWeights = cms.EDProducer("PdfWeightProducer",
        # Fix POWHEG if buggy (this PDF set will also appear on output,
        # so only two more PDF sets can be added in PdfSetNames if not "")
        #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
        GenTag = cms.untracked.InputTag(genParticleCollection),
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring(
            "CT10.LHgrid"
            #"cteq66.LHgrid"
            #, "MRST2006nnlo.LHgrid"
            #, "NNPDF10_100.LHgrid"
            #"cteq6mE.LHgrid"
            # ,"cteq6m.LHpdf"
            #"cteq6m.LHpdf"
        )
    )
else:
    process.pdfWeights = cms.Sequence()



####################################################################
## Event counter for events before any selection, including generator weights

from TopAnalysis.TopAnalyzer.CountEventAnalyzer_cfi import countEvents
process.EventsBeforeSelection = countEvents.clone()
process.EventsBeforeSelection.includePDFWeights = options.includePDFWeights
process.EventsBeforeSelection.pdfWeights = "pdfWeights:CT10"



####################################################################
## Write ntuple

from TopAnalysis.TopAnalyzer.NTupleWriter_cfi import writeNTuple

writeNTuple.sampleName = options.samplename
writeNTuple.channelName = options.mode
writeNTuple.systematicsName = options.systematicsName
writeNTuple.isMC = options.runOnMC
writeNTuple.isTtbarSample = topSignal
writeNTuple.isHiggsSample = higgsSignal
writeNTuple.isZSample = zGenInfo
writeNTuple.isMadgraphSample = madgraphSample

writeNTuple.includeTrigger = False
writeNTuple.includePdfWeights = options.includePDFWeights
writeNTuple.saveHadronMothers = False
writeNTuple.saveCHadronParticles = False

writeNTuple.electrons = isolatedElectronCollection
writeNTuple.muons = isolatedMuonCollection
writeNTuple.jets = jetCollection
writeNTuple.jetsForMet = jetForMetCollection
writeNTuple.jetsForMetUncorrected = jetForMetUncorrectedCollection
writeNTuple.met = metCollection
writeNTuple.mvaMet = mvaMetCollection
writeNTuple.vertices = selectedPrimaryVertices

writeNTuple.genParticles = genParticleCollection
writeNTuple.genJets = genJetCollection
writeNTuple.pdfWeights = "pdfWeights:CT10"

writeNTuple.BHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "BHadJetIndex")
writeNTuple.AntiBHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "AntiBHadJetIndex")
writeNTuple.BHadrons = cms.InputTag(genLevelBJetProducerInput, "BHadrons")
writeNTuple.AntiBHadrons = cms.InputTag(genLevelBJetProducerInput, "AntiBHadrons")
writeNTuple.BHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "BHadronFromTopB")
writeNTuple.AntiBHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronFromTopB")
writeNTuple.BHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "BHadronVsJet")
writeNTuple.AntiBHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronVsJet")

writeNTuple.genBHadPlusMothers = cms.InputTag(genBHadronMatcherInput, "genBHadPlusMothers")
writeNTuple.genBHadPlusMothersIndices = cms.InputTag(genBHadronMatcherInput, "genBHadPlusMothersIndices")
writeNTuple.genBHadIndex = cms.InputTag(genBHadronMatcherInput, "genBHadIndex")
writeNTuple.genBHadFlavour = cms.InputTag(genBHadronMatcherInput, "genBHadFlavour")
writeNTuple.genBHadJetIndex = cms.InputTag(genBHadronMatcherInput, "genBHadJetIndex")
writeNTuple.genBHadFromTopWeakDecay = cms.InputTag(genBHadronMatcherInput, "genBHadFromTopWeakDecay")
writeNTuple.genBHadLeptonIndex = cms.InputTag(genBHadronMatcherInput, "genBHadLeptonIndex")
writeNTuple.genBHadLeptonHadronIndex = cms.InputTag(genBHadronMatcherInput, "genBHadLeptonHadronIndex")
writeNTuple.genBHadLeptonViaTau = cms.InputTag(genBHadronMatcherInput, "genBHadLeptonViaTau")

writeNTuple.genCHadPlusMothers = cms.InputTag(genCHadronMatcherInput, "genCHadPlusMothers")
writeNTuple.genCHadPlusMothersIndices = cms.InputTag(genCHadronMatcherInput, "genCHadPlusMothersIndices")
writeNTuple.genCHadIndex = cms.InputTag(genCHadronMatcherInput, "genCHadIndex")
writeNTuple.genCHadFlavour = cms.InputTag(genCHadronMatcherInput, "genCHadFlavour")
writeNTuple.genCHadJetIndex = cms.InputTag(genCHadronMatcherInput, "genCHadJetIndex")
writeNTuple.genCHadFromTopWeakDecay = cms.InputTag(genCHadronMatcherInput, "genCHadFromTopWeakDecay")
writeNTuple.genCHadBHadronId = cms.InputTag(genCHadronMatcherInput, "genCHadBHadronId")
writeNTuple.genCHadLeptonIndex = cms.InputTag(genCHadronMatcherInput, "genCHadLeptonIndex")
writeNTuple.genCHadLeptonHadronIndex = cms.InputTag(genCHadronMatcherInput, "genCHadLeptonHadronIndex")
writeNTuple.genCHadLeptonViaTau = cms.InputTag(genCHadronMatcherInput, "genCHadLeptonViaTau")

process.writeNTuple = writeNTuple.clone()



####################################################################
## Path

process.path = cms.Path(
    process.pdfWeights *
    process.EventsBeforeSelection *
    process.topSplittingSequence *
    process.higgsGenSequence *
    process.zGenSequence *
    process.genJetSequence *
    process.genJetFlavourSequence *
    process.prefilterSequence *
    process.triggerSequence *
    process.goodOfflinePrimaryVertices *
    process.electronCorrectionSequence *
    process.correctRecoMuonEnergy *
    process.userPatSequence *
    process.mvaMetSequence *
    process.finalCollectionsSequence *
    process.preselectionSequence *
    process.jetProperties *
    process.madgraphWDecaySequence *
    process.writeNTuple
)

#pathnames = process.paths_().keys()



####################################################################
## Replace all input collections coherently to use corrections

# Corrected electrons
massSearchReplaceAnyInputTag(process.path,
                             cms.InputTag("gsfElectrons", ""),
                             cms.InputTag("calibratedElectrons", "calibratedGsfElectrons"),
                             True)
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons', '', 'RECO')
process.calibratedElectrons.inputElectronsTag = cms.InputTag('gsfElectrons', '', 'RECO')

# Corrected muons
massSearchReplaceAnyInputTag(process.path,
                             cms.InputTag("pfMuonsFromVertex"+pfpostfix),
                             cms.InputTag("correctMuonEnergy"),
                             True)
massSearchReplaceAnyInputTag(process.path,
                             cms.InputTag('muons'),
                             cms.InputTag("correctRecoMuonEnergy"),
                             True)
process.correctMuonEnergy.muonSrc = 'pfMuonsFromVertex'+pfpostfix
process.correctRecoMuonEnergy.muonSrc = 'muons'



####################################################################
## Signal catcher for more information on errors

process.load("TopAnalysis.TopUtils.SignalCatcher_cfi")



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
#process.path = cms.Path(process.printTree)
#process.pathNtuple = cms.Path()






