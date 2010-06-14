
import FWCore.ParameterSet.Config as cms

## ---
##   use this file to study different distributions for measurement of differential Cross Section
##   does also a monitoring of the cuts before they are applied
## ---

## ---
##    eventfilter is to get a special ttbar decay channel from ttbarSample by genmatching
##    decide whether to run on:
# 'background only' # 'all' # 'signal only' # 'semileptonic electron only' # 'dileptonic electron only' # 'dileptonic muon only' # 'fullhadronic' # 'dileptonic muon + electron only' # 'via single tau only' # 'dileptonic via tau only'
##    careful: genmatched selection- might cause problems for specific BG samples like qcd or data - use 'all' for them
##    signal is semileptonic with mu
##    background is ttbar other channels
##    'all' does no selection
## ---

eventFilter  = 'signal only'
## choose between # 'background only' # 'all' # 'signal only' # 'semileptonic electron only' # 'dileptonic electron only' # 'dileptonic muon only' # 'fullhadronic' # 'dileptonic muon + electron only' # 'via single tau only' # 'dileptonic via tau only'

writeOutput  = False # True

# analyse muon quantities
process = cms.Process("Selection")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(    

## add your favourite file here
    #'/store/user/henderle/Spring10/TTbar_MAD/PATtuple_10_1.root'
    #'/store/user/henderle/Spring10/TTbar_NLO/PATtuple_10_1.root'
    '/store/user/henderle/Spring10/WJets_MAD/PATtuple_100_2.root'
    #'/store/user/henderle/Spring10/ZJets_MAD/PATtuple_10_2.root'
    #'/store/user/henderle/Spring10/'
    )
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeDiffXSec_test.root')
)

# switch to run on data and remove all gen plots (type 'MC' or 'data')
if(not globals().has_key('runningOnData')): 
    runningOnData = "MC"

## ---
##    configure the cutflow scenario
## ---
## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
## high level trigger filter
process.load("TopAnalysis.TopFilter.sequences.triggerFilter_cff")
if(runningOnData == "MC"  ): process.hltMu9.TriggerResultsTag = "TriggerResults::REDIGI"
elif(runningOnData == "data"): print "" 
else: print "invalid choice of runningOnData"
## semileptonic selection
process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")
## generator matching
process.load("TopAnalysis.TopFilter.sequences.generatorMatching_cff")
## muon selection
process.load("TopAnalysis.TopFilter.sequences.muonSelection_cff")
## jet selection
process.load("TopAnalysis.TopFilter.sequences.jetSelection_cff")
## tool to select muons from gen Particles and save them as new collection
process.load("TopAnalysis.TopUtils.GenCandSelector_cfi")
## generator level based collections and semileptonic selection (muon and jets)
process.load("TopAnalysis.TopFilter.sequences.genSelection_cff")

## ---
## including analysis tools
## ---
## cross section module
process.load("TopAnalysis.TopAnalyzer.MuonCrossSection_cfi")
## jet kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.JetKinematics_cfi")
## muon kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.MuonKinematics_cfi")
## jet quality analyzer
process.load("TopAnalysis.TopAnalyzer.JetQuality_cfi")
## muon quality analyzer
process.load("TopAnalysis.TopAnalyzer.MuonQuality_cfi")
## MET analyzer
process.load("TopAnalysis.TopAnalyzer.METKinematics_cfi")
## electron kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.ElectronKinematics_cfi")
## muon jet kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.MuonJetKinematics_cfi")

## ---
##    set up vertex filter
## ---
process.PVSelection = cms.EDFilter("VertexSelector",
                                   src = cms.InputTag("offlinePrimaryVertices"),
                                   cut = cms.string("!isFake && ndof > 4 && abs(z) < 15 && position.Rho < 2"),
                                   filter = cms.bool(True),
                                   )

## ---
##    set up filter for different ttbar decay channels
## ---
process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
process.ttSemiLeptonicFilter = process.ttDecaySelection.clone()
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon = True
if(not eventFilter=='all'):
    ## adapt output filename
    if(eventFilter=='signal only'):
        process.TFileService.fileName = 'analyzeDiffXSec_testSig.root'
    elif(eventFilter=='background only'):
        process.ttSemiLeptonicFilter.invert = True
    elif(eventFilter=='semileptonic electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
    elif(eventFilter=='dileptonic electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
    elif(eventFilter=='dileptonic muon only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.muon     = True
    elif(eventFilter=='fullhadronic'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
    elif(eventFilter=='dileptonic muon + electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
    elif(eventFilter=='dileptonic via tau only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.tau      = True
    elif(eventFilter=='via single tau only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = True
    else:
        raise NameError, "'"+eventFilter+"' is not a prober eventFilter name choose: 'all', 'signal only', 'background only', 'semileptonic electron only', 'dileptonic electron only', 'dileptonic muon only', 'fullhadronic', 'via single tau only', 'dileptonic via tau only' or 'dileptonic muon + electron only'"
    
    ## sequence with filter for decay channel and trigger selection hltMu9
    process.filterSequence = cms.Sequence(process.makeGenEvt *
                                          process.ttSemiLeptonicFilter *
                                          process.hltMu9
                                          )
else:
    ## sequence without filter (only trigger selection hltMu9) - done when 'all' is chosen
    process.filterSequence = cms.Sequence(process.hltMu9)
    
## ---
##    set up genFilter for semileptonic muons and taus, where taus are decaying into leptons
## ---
process.ttSemiLeptonicFilterSemiTauMuon= process.ttDecaySelection.clone()
process.ttSemiLeptonicFilterSemiTauMuon.allowedTopDecays.decayBranchA.tau = True
process.ttSemiLeptonicFilterSemiTauMuon.allowedTopDecays.decayBranchA.muon= True
process.ttSemiLeptonicFilterSemiTauMuon.restrictTauDecays = cms.PSet(
    leptonic   = cms.bool(True),
    oneProng   = cms.bool(False),
    threeProng = cms.bool(False)
    )
process.genFilterSequence = cms.Sequence(  process.makeGenEvt
                                         + process.ttSemiLeptonicFilterSemiTauMuon)

## define ordered jets
uds0    = cms.PSet(index = cms.int32(0), correctionLevel = cms.string('abs') )
uds1    = cms.PSet(index = cms.int32(1), correctionLevel = cms.string('abs') )
uds2    = cms.PSet(index = cms.int32(2), correctionLevel = cms.string('abs') )
uds3    = cms.PSet(index = cms.int32(3), correctionLevel = cms.string('abs') )

## jet Kinematics to monitor JES shift
process.patLead_0_JetKinematics = process.analyzeJetKinematics.clone (src = 'selectedPatJets', analyze = uds0 )
process.patLead_1_JetKinematics = process.analyzeJetKinematics.clone (src = 'selectedPatJets', analyze = uds1 )
process.patLead_2_JetKinematics = process.analyzeJetKinematics.clone (src = 'selectedPatJets', analyze = uds2 )
process.patLead_3_JetKinematics = process.analyzeJetKinematics.clone (src = 'selectedPatJets', analyze = uds3 )
process.shiftedLead_0_JetKinematics = process.analyzeJetKinematics.clone (src = 'scaledJetEnergy:selectedPatJets', analyze = uds0 )
process.shiftedLead_1_JetKinematics = process.analyzeJetKinematics.clone (src = 'scaledJetEnergy:selectedPatJets', analyze = uds1 )
process.shiftedLead_2_JetKinematics = process.analyzeJetKinematics.clone (src = 'scaledJetEnergy:selectedPatJets', analyze = uds2 )
process.shiftedLead_3_JetKinematics = process.analyzeJetKinematics.clone (src = 'scaledJetEnergy:selectedPatJets', analyze = uds3 )

process.unshiftedJets = cms.Sequence(process.patLead_0_JetKinematics+
                                     process.patLead_1_JetKinematics+
                                     process.patLead_2_JetKinematics+
                                     process.patLead_3_JetKinematics
                                     )

process.shiftedJets = cms.Sequence(process.shiftedLead_0_JetKinematics+
                                   process.shiftedLead_1_JetKinematics+
                                   process.shiftedLead_2_JetKinematics+
                                   process.shiftedLead_3_JetKinematics
                                   )

## ---                                   
##    set up distribution for cross section measurement
## ---

## a) on reconstruction Niveau (with correlation plots for Migration effects if runningOnData=="MC")

if(runningOnData=="MC"):
    process.analyzeTightMuonCrossSectionRecNjets1 = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                      srcB = 'selectedGenMuonCollection')
    process.analyzeTightMuonCrossSectionRecNjets2 = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                      srcB = 'selectedGenMuonCollection')
    process.analyzeTightMuonCrossSectionRecNjets3 = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                      srcB = 'selectedGenMuonCollection')
    process.analyzeTightMuonCrossSectionRecNjets4 = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                      srcB = 'selectedGenMuonCollection')
    process.analyzeTightMuonCrossSectionRecNjets4Btag = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                          srcB = 'selectedGenMuonCollection')
    process.analyzeTightMuonCrossSectionRecNjets3Btag = process.analyzeCrossSectionMuonCorrelations.clone(srcA = 'tightMuons',
                                                                                                          srcB = 'selectedGenMuonCollection')

else:
    process.analyzeTightMuonCrossSectionRecNjets1     = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')
    process.analyzeTightMuonCrossSectionRecNjets2     = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')
    process.analyzeTightMuonCrossSectionRecNjets3     = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')
    process.analyzeTightMuonCrossSectionRecNjets4     = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')
    process.analyzeTightMuonCrossSectionRecNjets4Btag = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')
    process.analyzeTightMuonCrossSectionRecNjets3Btag = process.analyzeCrossSectionRecMuon.clone(srcA = 'tightMuons')   

## b) on generator Niveau
process.analyzeTightMuonCrossSectionGenNjets1 = process.analyzeCrossSectionGenMuon.clone(srcA = 'selectedGenMuonCollection')
process.analyzeTightMuonCrossSectionGenNjets2 = process.analyzeCrossSectionGenMuon.clone(srcA = 'selectedGenMuonCollection')
process.analyzeTightMuonCrossSectionGenNjets3 = process.analyzeCrossSectionGenMuon.clone(srcA = 'selectedGenMuonCollection')
process.analyzeTightMuonCrossSectionGenNjets4 = process.analyzeCrossSectionGenMuon.clone(srcA = 'selectedGenMuonCollection')

## ---
##    Set up selection steps for muon selection
## ---
process.combinedMuonsSelection        = process.muonSelection.clone (src = 'combinedMuons'       , minNumber = 1, maxNumber = 99999999)
process.triggerMuonsSelection         = process.muonSelection.clone (src = 'triggerMuons'        , minNumber = 1, maxNumber = 99999999)
process.trackMuonsSelection           = process.muonSelection.clone (src = 'trackMuons'          , minNumber = 1, maxNumber = 99999999)
process.goodMuonsSelection            = process.muonSelection.clone (src = 'goodMuons'           , minNumber = 1, maxNumber = 99999999)
process.goldenMuonsSelection          = process.muonSelection.clone (src = 'goldenMuons'         , minNumber = 1, maxNumber = 99999999)
process.tightMuonsSelection           = process.muonSelection.clone (src = 'tightMuons'          , minNumber = 1, maxNumber = 99999999)

process.muonCuts = cms.Sequence(process.combinedMuonsSelection        +
                                process.triggerMuonsSelection         +
                                process.trackMuonsSelection           +
                                process.goodMuonsSelection            +
                                process.goldenMuonsSelection          +
                                process.tightMuonsSelection           +
                                process.muonSelection
                                )

## ---
##    Set up selection steps for different jet multiplicities
## ---
process.leadingJetSelectionNjets1 = process.leadingJetSelection.clone (src = 'tightLeadingJets', minNumber = 1)
process.leadingJetSelectionNjets2 = process.leadingJetSelection.clone (src = 'tightLeadingJets', minNumber = 2)
process.leadingJetSelectionNjets3 = process.leadingJetSelection.clone (src = 'tightLeadingJets', minNumber = 3)
process.leadingJetSelectionNjets4 = process.leadingJetSelection.clone (src = 'tightLeadingJets', minNumber = 4)

## ---
##    collect selections for path 2 (jetmultiplicity 3 && btag) with different names
## ---
process.ttSemiLeptonicFilterb = process.ttSemiLeptonicFilter.clone()
process.hltMu9b = process.hltMu9.clone()
process.filterSequenceb = cms.Sequence(process.makeGenEvt *
                                       process.ttSemiLeptonicFilterb *
                                       process.hltMu9b
                                       )
process.bottomJetSelectionb        = process.bottomJetSelection.clone()
process.leadingJetSelectionNjets3b = process.leadingJetSelectionNjets3.clone()
process.muonSelectionb  = process.muonSelection.clone()
process.secondMuonVetob = process.secondMuonVeto.clone()
process.electronVetob   = process.electronVeto.clone()
process.PVSelectionb    = process.PVSelection.clone()

## ---
##    Set up selection steps for different (gen)-jet multiplicities
## ---
process.leadingGenJetSelectionNjets1 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 1)
process.leadingGenJetSelectionNjets2 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 2)
process.leadingGenJetSelectionNjets3 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 3)
process.leadingGenJetSelectionNjets4 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 4)

## ---
##    configure MET analyzer
## ---

process.analyzePfMET  = process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePatMET = process.analyzeMETCorrelations.clone(srcA = 'patMETs'  , srcB='tightMuons')
                                    
process.analyzePfMETNjets1    = process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePfMETNjets2    = process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePfMETNjets3    = process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePfMETNjets3Btag= process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePfMETNjets4    = process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')
process.analyzePfMETNjets4Btag= process.analyzeMETCorrelations.clone(srcA = 'patMETsPF', srcB='tightMuons')

process.analyzePatMETNjets1    = process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzePatMETNjets2    = process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzePatMETNjets3    = process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzePatMETNjets3Btag= process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzePatMETNjets4    = process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzePatMETNjets4Btag= process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')

## collect cut monitoring
## N-1 muon collections
process.noDbMuonQuality           = process.analyzeMuonQuality.clone   (src = 'noDbMuons'      )
process.noChi2MuonQuality         = process.analyzeMuonQuality.clone   (src = 'noChi2Muons'    )
process.noTrkHitsMuonQuality      = process.analyzeMuonQuality.clone   (src = 'noTrkHitsMuons' )
process.noIsoMuonQuality          = process.analyzeMuonQuality.clone   (src = 'noIsoMuons'     )
process.noEtaMuonKinematics       = process.analyzeMuonKinematics.clone(src = 'noEtaMuons'     )
process.noPtMuonKinematics        = process.analyzeMuonKinematics.clone(src = 'noPtMuons'      )
process.noDRMuonVetoJetsKinematics = process.analyzeMuonJetKinematics.clone(srcA = 'noDRMuons',
                                                                            srcB = 'vetoJets'  )
## N-1 jet collections
process.noEtaJetKinematics  = process.analyzeJetKinematics.clone(src = 'noEtaJets' )
process.noPtJetKinematics   = process.analyzeJetKinematics.clone(src = 'noPtJets'  )
process.noEmJetQuality      = process.analyzeJetQuality.clone(src = 'noEmJets'     )
process.noN90HitsJetQuality = process.analyzeJetQuality.clone(src = 'noN90HitsJets')
process.nofHPDJetQuality    = process.analyzeJetQuality.clone(src = 'nofHPDJets'   )
process.noPtLead_0_JetKinematics = process.analyzeJetKinematics.clone (src = 'noPtJets', analyze = uds0 )
process.noPtLead_1_JetKinematics = process.analyzeJetKinematics.clone (src = 'noPtJets', analyze = uds1 )
process.noPtLead_2_JetKinematics = process.analyzeJetKinematics.clone (src = 'noPtJets', analyze = uds2 )
process.noPtLead_3_JetKinematics = process.analyzeJetKinematics.clone (src = 'noPtJets', analyze = uds3 )
## veto collections
process.looseVetoMuonKinematics     = process.analyzeMuonKinematics.clone(src = 'looseMuons', analyze = cms.PSet(index = cms.int32(-1)))
process.looseVetoElectronKinematics = process.analyzeElectronKinematics.clone(src = 'looseElectrons',
                                                                              analyze = cms.PSet(index = cms.int32(-1)) )
process.patVetoElectronKinematics = process.analyzeElectronKinematics.clone(src = 'selectedPatElectrons',
                                                                            analyze = cms.PSet(index = cms.int32(-1)) )
## muon cutflow
process.combinedMuonKinematics = process.analyzeMuonKinematics.clone(src = 'combinedMuons')
process.triggerMuonQuality     = process.analyzeMuonQuality.clone   (src = 'triggerMuons' )
process.trackMuonKinematics    = process.analyzeMuonKinematics.clone(src = 'trackMuons'   )
process.goodMuonVetoJetsKinematics = process.analyzeMuonJetKinematics.clone(srcA = 'goodMuons',
                                                                            srcB = 'vetoJets'  )
process.goldenMuonQuality      = process.analyzeMuonQuality.clone   (src = 'goldenMuons'  )
process.tightMuonKinematics    = process.analyzeMuonKinematics.clone(src = 'tightMuons'   )
process.tightMuonQuality       = process.analyzeMuonQuality.clone   (src = 'tightMuons'   )

## jet cutflow
process.patJetKinematics = process.analyzeJetKinematics.clone(src = 'selectedPatJets')
process.centralLead_0_JetKinematics = process.analyzeJetKinematics.clone (src = 'centralJets', analyze = uds0 )
process.centralLead_1_JetKinematics = process.analyzeJetKinematics.clone (src = 'centralJets', analyze = uds1 )
process.centralLead_2_JetKinematics = process.analyzeJetKinematics.clone (src = 'centralJets', analyze = uds2 )
process.centralLead_3_JetKinematics = process.analyzeJetKinematics.clone (src = 'centralJets', analyze = uds3 )
process.centralJetKinematics = process.analyzeJetKinematics.clone(src = 'centralJets'    )
process.reliableJetQuality   = process.analyzeJetQuality.clone   (src = 'reliableJets'   )
process.tightJetKinematics  = process.analyzeJetKinematics.clone(src = 'tightLeadingJets')
## btag selection cuts
process.tightJetQuality     = process.analyzeJetQuality.clone   (src = 'tightLeadingJets')
process.bottomJetKinematics = process.analyzeJetKinematics.clone(src = 'tightBottomJets' )


process.monitorNMinusOneMuonCuts = cms.Sequence(process.noDbMuonQuality            +
                                                process.noChi2MuonQuality          +
                                                process.noTrkHitsMuonQuality       +
                                                process.noIsoMuonQuality           +
                                                process.noEtaMuonKinematics        +
                                                process.noPtMuonKinematics         +
                                                process.noDRMuonVetoJetsKinematics
                                                )
process.monitorMuonCutflow = cms.Sequence(process.combinedMuonKinematics     +
                                          process.triggerMuonQuality         +
                                          process.trackMuonKinematics        +
                                          process.goodMuonVetoJetsKinematics +
                                          process.goldenMuonQuality          +
                                          process.tightMuonKinematics        +
                                          process.tightMuonQuality           +
                                          process.analyzePatMET              +
                                          process.analyzePfMET               
                                          )
process.monitorNMinusOneJetCuts = cms.Sequence(process.noPtLead_0_JetKinematics  +
                                               process.noPtLead_1_JetKinematics  +
                                               process.noPtLead_2_JetKinematics  +
                                               process.noPtLead_3_JetKinematics  +
                                               process.noEtaJetKinematics        +
                                               process.noPtJetKinematics         +
                                               process.noEmJetQuality            +
                                               process.noN90HitsJetQuality       +
                                               process.nofHPDJetQuality       
                                               )
process.monitorJetCutflow = cms.Sequence(process.patJetKinematics            +
                                         process.centralLead_0_JetKinematics +
                                         process.centralLead_1_JetKinematics +
                                         process.centralLead_2_JetKinematics +
                                         process.centralLead_3_JetKinematics +
                                         process.centralJetKinematics        +
                                         process.reliableJetQuality          +
                                         process.tightJetKinematics
                                         )
process.monitorVetoCuts = cms.Sequence(process.looseVetoMuonKinematics     +
                                       process.patVetoElectronKinematics   +
                                       process.looseVetoElectronKinematics
                                       )
process.monitorBtagCuts = cms.Sequence(process.tightJetQuality     +
                                       process.bottomJetKinematics 
                                       )

## ---
##     collect analyzers for different jet multiplicities
##    ( Selection - CrossSectionMeasurement - CorrelationPlots - JetKinematics - METKinematics)
## ---
process.jetMultiplicity1 = cms.Sequence(process.leadingJetSelectionNjets1             +
                                        process.analyzeTightMuonCrossSectionRecNjets1 +
                                        process.analyzePfMETNjets1                    +
                                        process.analyzePatMETNjets1                   )
process.jetMultiplicity2 = cms.Sequence(process.leadingJetSelectionNjets2             +
                                        process.analyzeTightMuonCrossSectionRecNjets2 +
                                        process.analyzePfMETNjets2                    +
                                        process.analyzePatMETNjets2                   )
process.jetMultiplicity3 = cms.Sequence(process.leadingJetSelectionNjets3             +
                                        process.analyzeTightMuonCrossSectionRecNjets3 +
                                        process.analyzePfMETNjets3                    +
                                        process.analyzePatMETNjets3                   )
process.jetMultiplicity4 = cms.Sequence(process.leadingJetSelectionNjets4             +
                                        process.analyzeTightMuonCrossSectionRecNjets4 +
                                        process.analyzePfMETNjets4                    +
                                        process.analyzePatMETNjets4                   )
process.jetMultiplicity3Btag = cms.Sequence(process.leadingJetSelectionNjets3b                 +
                                            process.bottomJetSelectionb                        +
                                            process.analyzeTightMuonCrossSectionRecNjets3Btag  +
                                            process.analyzePfMETNjets3Btag                     +
                                            process.analyzePatMETNjets3Btag                    )
process.jetMultiplicity4Btag = cms.Sequence(process.bottomJetSelection                         +
                                            process.analyzeTightMuonCrossSectionRecNjets4Btag  +
                                            process.analyzePfMETNjets4Btag                     +
                                            process.analyzePatMETNjets4Btag                    )

# ## produce decaySubset
# process.load("TopQuarkAnalysis.TopEventProducers.producers.TopDecaySubset_cfi") 
# ## add message logger
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# #process.MessageLogger.categories.append(*'TopDecaySubset_printTarget'*)
# process.MessageLogger.categories.append('TopDecaySubset_printSource')
# process.MessageLogger.cerr.TopDecaySubset_printSource = cms.untracked.PSet(
#    limit = cms.untracked.int32(10)
# )


## ---
##    run the final sequences
## ---

process.p1 = cms.Path(
                      ## do the gen event selection (decay channel) and the trigger selection (hltMu9)
                      process.filterSequence                        *
                      ## do the PV event selection
                      process.PVSelection                           *
                      ## introduce some collections
                      process.semiLeptonicSelection                 *
                      process.selectNMinusOneJets                   *
                      process.selectNMinusOneMuons                  *
                      process.isolatedGenMuons                      *
                      process.semiLeptGenCollections                *
                      ## monitor all muon cut quantities
                      process.monitorNMinusOneMuonCuts              *
                      process.monitorMuonCutflow                    *
                      ## do the event selection for muon
                      process.muonCuts                              *
                      ## monitor veto collection
                      process.monitorVetoCuts                       *
                      ## do event selection veto cuts
                      process.secondMuonVeto                        *
                      process.electronVeto                          *
                      ## monitor all jet cut quantities
                      process.monitorNMinusOneJetCuts               *
                      process.monitorJetCutflow                     *
                      ## analysis for different jet multiplicities
                      ## N_jets >= 1                    
                      process.jetMultiplicity1                      *
                      ## N_jets >= 2
                      process.jetMultiplicity2                      *
                      ## N_jets >= 3
                      process.jetMultiplicity3                      *
                      ## N_jets >= 4
                      process.jetMultiplicity4                      *
                      ## monitor btag cut quantities
                      process.monitorBtagCuts                       *
                      ## N_jets >= 4 + btag >=1                      
                      process.jetMultiplicity4Btag
                      )
## Njets>=3 & btag
process.p2 = cms.Path(
                      ## do the gen event selection (decay channel) and the trigger selection (hltMu9)
                      process.filterSequenceb                       *
                      ## do the PV event selection
                      process.PVSelectionb                          *
                      ## introduce some collections
                      process.semiLeptonicSelection                 *
                      process.isolatedGenMuons                      *
                      process.semiLeptGenCollections                *
                      ## do the event selection for muon
                      process.muonSelectionb                        *
                      ## do event selection veto cuts
                      process.secondMuonVetob                       *
                      process.electronVetob                         *
                      ## N_jets >= 3 + btag >=1
                      process.jetMultiplicity3Btag
                      )
## on generator niveau
if(runningOnData=="MC"):
    print "running on Monte Carlo, gen-plots produced"
    process.p3 = cms.Path(
        ## gen event selection: semileptonic (muon & tau->lepton)
        process.genFilterSequence                     *
        ## introduce some collections
        process.isolatedGenMuons                      *
        process.semiLeptGenCollections                *
        ## do the event selection for muon
        process.genMuonSelection                      *
     #   process.decaySubset                           *
        ## for N_jets = 1+
        process.leadingGenJetSelectionNjets1          *
        process.analyzeTightMuonCrossSectionGenNjets1 *
        ## for N_jets = 2+
        process.leadingGenJetSelectionNjets2          *
        process.analyzeTightMuonCrossSectionGenNjets2 *
        ## for N_jets = 3+
        process.leadingGenJetSelectionNjets3          *
        process.analyzeTightMuonCrossSectionGenNjets3 *
        ##  for N_jets = 4+
        process.leadingGenJetSelectionNjets4          *
        process.analyzeTightMuonCrossSectionGenNjets4 
        )
elif(runningOnData=="data"):
    print "running on data, no gen-plots"
else:
    print "choose runningOnData= data or MC, creating no gen-plots"

## Output Module Configuration
if(writeOutput):
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('patTuple_selected.root'),
                                   # save only events passing the full path
                                   SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p1') ),
                                   #save output (comment to keep everything...)
                                   #outputCommands = cms.untracked.vstring('drop *',) 
                                   )
    process.outpath = cms.EndPath(process.out)
    
