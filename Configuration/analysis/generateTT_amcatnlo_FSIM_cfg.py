import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
import glob

#inputfiles = []
#hepmcfiles = glob.glob('*.lhe')
#for hepmcfile in hepmcfiles:
#  inputfiles.append('file:' + hepmcfile)
#hepmcfiles = glob.glob('../../prolog/*.lhe')
#for hepmcfile in hepmcfiles:
#  inputfiles.append('file:' + hepmcfile)
#print inputfiles

options = VarParsing.VarParsing ('standard')

options.register('frag', 'pythia8', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "parton shower and hadronization (pythia8|herwigpp)")
options.register('myjobid', 0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "job id")
options.register('maxevents', 10000, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "events per job")


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

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_2012_Summer_inTimeOnly_cff')
process.load('FastSimulation.Configuration.Geometries_START_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('HLTrigger.Configuration.HLT_7E33v2_Famos_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxevents)
)
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# Input source
#process.source = cms.Source("LHESource",
#    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring(inputfiles)
#)
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

process.source.firstEvent =  cms.untracked.uint32(options.myjobid * options.maxevents + 1)
process.source.firstLuminosityBlock =  cms.untracked.uint32(options.myjobid + 1)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.simulation = cms.Sequence(process.simulationWithFamos)
process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)
process.Realistic8TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V7C::All', '')

process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    nEvents = cms.uint32(options.maxevents),
    args = cms.vstring('slc6_amd64_gcc472/8TeV/madgraph/V5_2.2.1/tt_5f_NLO_'+options.frag+'/v1/', 
        'tt_5f_NLO_'+options.frag),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball.sh'),
    numberOfParameters = cms.uint32(2),
    outputFile = cms.string('tt_5f_NLO_'+options.frag+'_final.lhe')
)


if (options.frag == 'pythia8'):
    process.generator = cms.EDFilter("Pythia8HadronizerFilter",
        maxEventsToPrint = cms.untracked.int32(1),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        filterEfficiency = cms.untracked.double(1.0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        comEnergy = cms.double(8000.),
        PythiaParameters = cms.PSet(
        processParameters = cms.vstring(
            'Main:timesAllowErrors = 10000',
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tau0Max = 10',
            'Tune:ee 3',
            'Tune:pp 5',
            'SpaceShower:pTmaxMatch = 1',
            'SpaceShower:pTmaxFudge = 1',
            'SpaceShower:MEcorrections = off',
            'TimeShower:pTmaxMatch = 1',
            'TimeShower:pTmaxFudge = 1',
            'TimeShower:MEcorrections = off',
            'TimeShower:globalRecoil = on',
            'TimeShower:limitPTmaxGlobal = on',
            'TimeShower:nMaxGlobalRecoil = 1',
            'TimeShower:globalRecoilMode = 2',
            'TimeShower:nMaxGlobalBranch = 1',
            'SLHA:keepSM = on',
            'SLHA:minMassSM = 1000.',
            'Check:epTolErr = 0.01',
            '6:m0 = 172.5',    # top mass'
            ),
            parameterSets = cms.vstring('processParameters')
        )
    )
elif (options.frag == 'herwigpp'):
    from Configuration.Generator.HerwigppDefaults_cfi import *
    from Configuration.Generator.HerwigppUE_EE_3C_cfi import *
    process.generator = cms.EDFilter("LHEProducer",
        eventsToPrint = cms.untracked.uint32(1),

        hadronisation = cms.PSet(
            herwigDefaultsBlock,
            herwigppUESettingsBlock,
            #herwigValidationBlock,

            generator = cms.string('ThePEG'),

            configFiles = cms.vstring(),

            parameterSets = cms.vstring(
                'basicSetup',
                'herwigppUE_EE_3C_8000GeV',
                
                #'setParticlesStableForDetector',
                'lheDefaults', 
                'productionParameters',
                #'lheDefaultPDFs'
                #'powhegDefaults',
            ),
            productionParameters = cms.vstring(
                # 1.) NECESSARY SETTINGS FOR RUNNING WITH MC@NLO EVENTS (DO NOT MODIFY)
                'set /Herwig/Shower/Evolver:HardVetoMode 1',
                'set /Herwig/Shower/Evolver:HardVetoScaleSource 1',
                'set /Herwig/Shower/Evolver:MECorrMode 0',
                
                #Boost and reconstruction stuff
                'set /Herwig/Shower/KinematicsReconstructor:ReconstructionOption General',
                'set /Herwig/Shower/KinematicsReconstructor:InitialInitialBoostOption LongTransBoost',

                # create the Handler & Reader
                'set /Herwig/EventHandlers/LHEReader:AllowedToReOpen 0',
                'set /Herwig/EventHandlers/LHEReader:MomentumTreatment RescaleEnergy',
                'set /Herwig/EventHandlers/LHEReader:WeightWarnings 0',

                'set /Herwig/EventHandlers/LHEHandler:WeightOption VarNegWeight',
                'set /Herwig/EventHandlers/LHEHandler:PartonExtractor /Herwig/Partons/QCDExtractor',
                'set /Herwig/EventHandlers/LHEHandler:CascadeHandler /Herwig/Shower/ShowerHandler ',
                'set /Herwig/EventHandlers/LHEHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
                'set /Herwig/EventHandlers/LHEHandler:DecayHandler /Herwig/Decays/DecayHandler',
                'insert /Herwig/EventHandlers/LHEHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler',

                'set /Herwig/Generators/LHCGenerator:PrintEvent 1',
                'set /Herwig/Generators/LHCGenerator:DebugLevel 1',
                
                # Use internal HW++ defualt PDF sets
                'set /Herwig/EventHandlers/LHEReader:PDFA /Herwig/Partons/MRST-NLO',
                'set /Herwig/EventHandlers/LHEReader:PDFB /Herwig/Partons/MRST-NLO',
                'set /Herwig/Particles/p+:PDF /Herwig/Partons/MRST-NLO',
                'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/MRST-NLO',
                'set /Herwig/Shower/ShowerHandler:PDFA /Herwig/Partons/MRST-NLO',
                'set /Herwig/Shower/ShowerHandler:PDFB /Herwig/Partons/MRST-NLO',
                
                # top mass
                'set /Herwig/Particles/t:NominalMass 172.5*GeV',
                'set /Herwig/Particles/tbar:NominalMass 172.5*GeV',
                
                'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV',
            ),
        )
    )

# Output definition
commands = process.AODSIMEventContent.outputCommands

commands.extend([
#'keep *',
'drop CaloTowersSorted_towerMaker__HLT',
#'drop recoPFCandidates_FSparticleFlow__HLT',
#'drop EcalRecHitsSorted_reducedEcalRecHitsEE__HLT',
'drop triggerTriggerEvent_hltTriggerSummaryAOD__HLT',
'drop recoPFTaus_hpsPFTauProducer__HLT',
#'drop recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_HLT',
'drop recoPFJets_kt6PFJets__HLT',
'drop recoPFJets_kt4PFJets__HLT',
'drop recoPFJets_ak5PFJets__HLT',
'drop recoPFJets_ak7PFJets__HLT',
'drop recoTrackExtrapolations_trackExtrapolator__HLT',
#'drop EcalRecHitsSorted_reducedEcalRecHitsEB__HLT',
'drop recoCaloJets_kt6CaloJets__HLT',
#'drop recoCaloJets_ak5CaloJets__HLT',
'drop recoCaloJets_ak7CaloJets__HLT',
'drop recoCaloJets_kt4CaloJets__HLT',
#'drop recoMuons_muons__HLT',
'drop recoJPTJets_JetPlusTrackZSPCorJetAntiKt5__HLT',
#'drop recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower__HLT',
'drop recoPhotons_photons__HLT',
'drop recoJetIDedmValueMap_kt6JetID__HLT',
'drop recoJetIDedmValueMap_ak5JetID__HLT',
'drop recoJetIDedmValueMap_kt4JetID__HLT',
'drop recoJetIDedmValueMap_ak7JetID__HLT',
'drop recoPhotons_pfPhotonTranslator_pfphot_HLT',
'drop recoConversions_allConversions__HLT',
#'drop recoCaloClusters_hybridSuperClusters_hybridBarrelBasicClusters_HLT',
#'drop recoSuperClusters_correctedHybridSuperClusters__HLT',
'drop recoJetedmRefToBaseProdTorecoJetExtendedAssociationJetExtendedDatasAssociationVector_kt4JetExtender__HLT',
'drop recoJetedmRefToBaseProdTorecoJetExtendedAssociationJetExtendedDatasAssociationVector_ak5JetExtender__HLT',
#'drop recoPFCandidates_particleFlowTmp_electrons_HLT',
'drop recoJetedmRefToBaseProdTorecoJetExtendedAssociationJetExtendedDatasAssociationVector_ak7JetExtender__HLT',
'drop recoGenJets_ak7GenJets__HLT',
'drop recoGenJets_kt6GenJets__HLT',
#'drop recoGsfTracks_electronGsfTracks__HLT',
'drop recoGenJets_kt4GenJets__HLT',
'drop TrackingRecHitsOwned_refittedStandAloneMuons__HLT',
'drop TrackingRecHitsOwned_standAloneMuons__HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedPunchThroughNeutralHadrons_HLT',
#'drop recoIsoDepositedmValueMap_muIsoDepositCalByAssociatorTowers_ecal_HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedTrackerAndGlobalMuons_HLT',
'drop recoVertexCompositeCandidates_generalV0Candidates_Kshort_HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedPunchThroughMuons_HLT',
#'drop recoPFCandidates_particleFlowTmp_AddedMuonsAndHadrons_HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedCosmicsMuons_HLT',
'drop recoJetedmRefToBaseProdrecoTracksrecoTrackrecoTracksTorecoTrackedmrefhelperFindUsingAdvanceedmRefVectorsAssociationVector_ak5JetTracksAssociatorAtVertex__HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedFakeMuons_HLT',
#'drop recoIsoDepositedmValueMap_muons_ecal_HLT',
'drop recoCaloMETs_corMetGlobalMuons__HLT',
#'drop recoPFCandidates_particleFlowTmp_CleanedHF_HLT',
'drop recoCaloMETs_metHO__HLT',
'drop recoCaloMETs_met__HLT',
'drop recoCaloMETs_metNoHF__HLT',
'drop recoJetedmRefToBaseProdrecoTracksrecoTrackrecoTracksTorecoTrackedmrefhelperFindUsingAdvanceedmRefVectorsAssociationVector_ak7JetTracksAssociatorAtVertex__HLT',
'drop recoTrackJets_ak5TrackJets__HLT',
'drop recoPreshowerClusters_multi5x5SuperClustersWithPreshower_preshowerYClusters_HLT',
#'drop recoIsoDepositedmValueMap_muons_muIsoDepositTk_HLT',
#'drop recoIsoDepositedmValueMap_muIsoDepositCalByAssociatorTowers_hcal_HLT',
#'drop recoIsoDepositedmValueMap_muIsoDepositTk__HLT',
'drop recoTrackExtras_refittedStandAloneMuons__HLT',
'drop recoRecoTauPiZeros_hpsPFTauProducer_pizeros_HLT',
'drop recoGenMETs_genMetCaloAndNonPrompt__HLT',
'drop recoMuonTimeExtraedmValueMap_muons_combined_HLT',
'drop recoConversions_pfPhotonTranslator_pfphot_HLT',
#'drop recoGenMETs_genMetTrue__HLT',
'drop recoGenMETs_genMetCalo__HLT',
#'drop recoIsoDepositedmValueMap_muons_hcal_HLT',
'drop recoMuonTimeExtraedmValueMap_muons_dt_HLT',
'drop recoMuonTimeExtraedmValueMap_muons_csc_HLT',
#'drop recoIsoDepositedmValueMap_muIsoDepositCalByAssociatorTowers_ho_HLT',
#'drop recoIsoDepositedmValueMap_muons_muIsoDepositJets_HLT',
'drop recoVertexCompositeCandidates_generalV0Candidates_Lambda_HLT',
#'drop recoIsoDepositedmValueMap_muIsoDepositJets__HLT',
#'drop recoPreshowerClusters_multi5x5SuperClustersWithPreshower_preshowerXClusters_HLT',
#'drop recoSuperClusters_pfElectronTranslator_pf_HLT',
'drop recoRecoChargedRefCandidates_trackRefsForJets__HLT',
#'drop EcalRecHitsSorted_reducedEcalRecHitsES__HLT',
#'drop recoIsoDepositedmValueMap_muons_ho_HLT',
'drop l1extraL1MuonParticles_l1extraParticles__HLT',
#'drop recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_HLT',
#'drop recoBeamSpot_offlineBeamSpot__HLT',
#'drop recoTracks_refittedStandAloneMuons_UpdatedAtVtx_HLT',
#'drop recoTracks_standAloneMuons_UpdatedAtVtx_HLT',
#'drop recoTracks_refittedStandAloneMuons__HLT',
'drop recoMuonShoweredmValueMap_muons_muonShowerInformation_HLT',
#'drop recoSuperClusters_pfPhotonTranslator_pfphot_HLT',
#'drop recoTracks_tevMuons_default_HLT',
#'drop recoTracks_tevMuons_picky_HLT',
#'drop recoTracks_tevMuons_firstHit_HLT',
#'drop recoTracks_globalMuons__HLT',
#'drop recoTracks_standAloneMuons__HLT',
'drop l1extraL1EtMissParticles_l1extraParticles_MET_HLT',
'drop l1extraL1EtMissParticles_l1extraParticles_MHT_HLT',
#'drop recoTracks_tevMuons_dyt_HLT',
'drop recoMuonMETCorrectionDataedmValueMap_muonTCMETValueMapProducer_muCorrData_HLT',
'drop recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_muCorrData_HLT',
#'drop recoSuperClusters_hfEMClusters__HLT',
'drop LHEEventProduct_externalLHEProducer__LHE',
'drop recoMETs_tcMetWithPFclusters__HLT',
#'drop recoPFRecHits_particleFlowClusterHFHAD_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowClusterECAL_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowClusterHFEM_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowClusterHCAL_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowRecHitECAL_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowRecHitHCAL_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowClusterHO_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowClusterPS_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowRecHitHO_Cleaned_HLT',
#'drop recoPFRecHits_particleFlowRecHitPS_Cleaned_HLT',
'drop recoMETs_htMetAK7__HLT',
'drop recoMETs_htMetKT6__HLT',
'drop recoMETs_htMetKT4__HLT',
'drop recoMETs_htMetIC5__HLT',
'drop recoMETs_htMetAK5__HLT',
'drop recoMETs_tcMet__HLT',
#'drop l1extraL1EmParticles_l1extraParticles_NonIsolated_HLT',
#'drop recoPFCandidateedmPtredmValueMap_particleFlow_electrons_HLT',
'drop l1extraL1JetParticles_l1extraParticles_Central_HLT',
#'drop recoPFCandidateedmPtredmValueMap_particleFlow_muons_HLT',
'drop l1extraL1JetParticles_l1extraParticles_Forward_HLT',
#'drop l1extraL1EmParticles_l1extraParticles_Isolated_HLT',
#'drop recoPFCandidateedmPtredmValueMap_particleFlow_photons_HLT',
'drop l1extraL1JetParticles_l1extraParticles_Tau_HLT',
#'drop recoCaloClusters_pfElectronTranslator_pf_HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3rawElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2rawElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByIsolationMVA2raw__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByIsolationMVAraw__HLT',
#'drop GenEventInfoProduct_generator__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3rawElectronRejection_category_HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2rawElectronRejection_category_HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3VTightElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3MediumElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2VLooseElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2MediumElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3TightElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2TightElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA2LooseElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVA3LooseElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightMuonRejection2__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightMuonRejection3__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseMuonRejection3__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByDeadECALElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightMuonRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMVAElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumChargedIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseElectronRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByVLooseChargedIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightChargedIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseChargedIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumMuonRejection2__HLT',
#'drop recoSuperClustersToOnerecoHFEMClusterShapesAssociation_hfEMClusters__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumIsolationMVA2__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseMuonRejection2__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumMuonRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightIsolationMVA2__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumIsolationMVA__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseIsolationMVA2__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseMuonRejection__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightIsolationMVA__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseIsolationMVA__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByDecayModeFinding__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByVLooseIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByMediumIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByLooseIsolation__HLT',
'drop recoPFTauDiscriminator_hpsPFTauDiscriminationByTightIsolation__HLT',
'drop recoTracksToOnerecoTracksAssociation_tevMuons_firstHit_HLT',
'drop recoTracksToOnerecoTracksAssociation_tevMuons_default_HLT',
'drop HBHERecHitsSorted_reducedHcalRecHits_hbhereco_HLT',
'drop recoTracksToOnerecoTracksAssociation_tevMuons_picky_HLT',
'drop recoTracksToOnerecoTracksAssociation_tevMuons_dyt_HLT',
'drop recoPhotonCores_photonCore__HLT',
'drop recoRecoEcalCandidates_hfRecoEcalCandidate__HLT',
'drop recoPhotonCores_pfPhotonTranslator_pfphot_HLT',
'drop HORecHitsSorted_reducedHcalRecHits_horeco_HLT',
'drop HFRecHitsSorted_reducedHcalRecHits_hfreco_HLT',
#'drop recoCaloClusters_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters_HLT',
'drop recoHFEMClusterShapes_hfEMClusters__HLT',
'drop recoJetedmRefToBaseProdTofloatsAssociationVector_combinedSecondaryVertexMVABJetTags__HLT',
#'drop recoCaloClusters_pfPhotonTranslator_pfphot_HLT',
#'drop recoCaloClusters_hfEMClusters__HLT',
'drop L1GlobalTriggerReadoutRecord_gtDigis__HLT',
'drop l1extraL1HFRingss_l1extraParticles__HLT',
])

commands.extend([
#'drop recoPreshowerClusters_pfElectronTranslator_pf_HLT',
'drop recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_multi5x5PreshowerYClustersShape_HLT',
'drop recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_multi5x5PreshowerXClustersShape_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueChargedAll04_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueChargedAll03_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueCharged04_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueCharged03_HLT',
#'drop recoPreshowerClusters_pfPhotonTranslator_pfphot_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueGamma04_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValuePU04_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueGamma03_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValuePU03_HLT',
#'drop booledmValueMap_PhotonIDProd_PhotonCutBasedIDLooseEM_HLT',
#'drop booledmValueMap_PhotonIDProd_PhotonCutBasedIDLoose_HLT',
#'drop booledmValueMap_PhotonIDProd_PhotonCutBasedIDTight_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueNeutral04_HLT',
#'drop doubleedmValueMap_muons_muPFIsoValueNeutral03_HLT',
'drop double_kt6PFJetsCentralChargedPileUp_sigma_HLT',
'drop double_kt6PFJetsCentralChargedPileUp_rho_HLT',
'drop floatedmValueMap_eidRobustHighEnergy__HLT',
'drop double_kt6PFJetsCentralNeutral_sigma_HLT',
'drop floatedmValueMap_eidRobustLoose__HLT',
'drop floatedmValueMap_eidRobustTight__HLT',
'drop double_kt6PFJetsCentralNeutral_rho_HLT',
'drop double_kt6CaloJetsCentral_sigma_HLT',
'drop floatedmValueMap_eidLoose__HLT',
'drop double_fixedGridRhoFastjetAll__HLT',
'drop floatedmValueMap_eidTight__HLT',
'drop double_kt6PFJetsCentralNeutralTight_sigma_HLT',
'drop double_kt6CaloJetsCentral_rho_HLT',
'drop double_kt6PFJetsCentralNeutralTight_rho_HLT',
'drop double_kt6CaloJets_sigma_HLT',
'drop double_kt6CaloJets_rho_HLT',
'drop double_kt6PFJets_sigma_HLT',
'drop double_kt6PFJets_rho_HLT',
'drop doubles_ak5TrackJets_sigmas_HLT',
'drop doubles_kt4CaloJets_sigmas_HLT',
'drop doubles_kt6CaloJets_sigmas_HLT',
'drop doubles_ak5CaloJets_sigmas_HLT',
'drop doubles_ak7CaloJets_sigmas_HLT',
'drop double_fixedGridRhoAll__HLT',
'drop doubles_kt4GenJets_sigmas_HLT',
'drop doubles_kt6GenJets_sigmas_HLT',
'drop doubles_ak7GenJets_sigmas_HLT',
'drop doubles_ak5GenJets_sigmas_HLT',
'drop doubles_ak5TrackJets_rhos_HLT',
'drop doubles_ak5PFJets_sigmas_HLT',
'drop doubles_kt4PFJets_sigmas_HLT',
'drop doubles_kt6PFJets_sigmas_HLT',
'drop doubles_ak7PFJets_sigmas_HLT',
'drop doubles_ak5CaloJets_rhos_HLT',
'drop doubles_kt6CaloJets_rhos_HLT',
'drop doubles_ak7CaloJets_rhos_HLT',
'drop doubles_kt4CaloJets_rhos_HLT',
'drop doubles_kt6GenJets_rhos_HLT',
'drop doubles_ak5GenJets_rhos_HLT',
'drop doubles_ak7GenJets_rhos_HLT',
'drop doubles_kt4GenJets_rhos_HLT',
'drop doubles_ak5PFJets_rhos_HLT',
'drop doubles_ak7PFJets_rhos_HLT',
'drop doubles_kt4PFJets_rhos_HLT',
'drop double_ak5TrackJets_sigma_HLT',
'drop double_ak7CaloJets_sigma_HLT',
'drop double_kt4CaloJets_sigma_HLT',
'drop double_ak5CaloJets_sigma_HLT',
'drop double_ak5GenJets_sigma_HLT',
'drop double_ak7GenJets_sigma_HLT',
'drop double_kt6GenJets_sigma_HLT',
'drop double_kt4GenJets_sigma_HLT',
'drop double_ak5TrackJets_rho_HLT',
'drop double_ak5PFJets_sigma_HLT',
'drop double_ak7PFJets_sigma_HLT',
'drop double_kt4PFJets_sigma_HLT',
'drop double_kt4CaloJets_rho_HLT',
'drop double_ak5CaloJets_rho_HLT',
'drop double_ak7CaloJets_rho_HLT',
'drop double_ak7GenJets_rho_HLT',
'drop double_ak5GenJets_rho_HLT',
'drop double_kt4GenJets_rho_HLT',
'drop double_kt6GenJets_rho_HLT',
'drop double_ak5PFJets_rho_HLT',
'drop double_ak7PFJets_rho_HLT',
'drop double_kt4PFJets_rho_HLT',
])

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = commands, #process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:TT_FSIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AODSIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Path and EndPath definitions
process.lhe_step = cms.Path(process.externalLHEProducer)
process.generation_step = cms.Path(process.pgen_genonly)
process.reconstruction = cms.Path(process.reconstructionWithFamos)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.lhe_step,process.generation_step,process.genfiltersummary_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.reconstruction,process.AODSIMoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
        if path in ['lhe_step']: continue
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
