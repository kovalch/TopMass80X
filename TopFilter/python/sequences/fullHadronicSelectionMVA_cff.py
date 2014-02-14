import FWCore.ParameterSet.Config as cms

## all leptons kept for possible rejection of events later on (not yet implemented)

## jet selector
from TopAnalysis.TopFilter.sequences.jetSelection_cff import *
## muon selector
#from TopAnalysis.TopFilter.sequences.muonSelection_cff import *

## jet selector
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
## muon selector
#from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
## electron selector
#from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *

## jet count filter
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
## muon count filter
#from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
## electron count filter
#from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *

## ---
##    setup the collections for the fully-hadronic event selection
## ---

## define tight cuts for jet kinematics
#tight4JetCut  = 'pt > 60. & abs(eta) < 2.4'
tight4JetCut  = 'pt > 30. & abs(eta) < 2.4'
tight5JetCut  = 'pt > 30. & abs(eta) < 2.4'
tight6JetCut  = 'pt > 30. & abs(eta) < 2.4'
#tight5JetCut  = 'pt > 50. & abs(eta) < 2.4'
#tight6JetCut  = 'pt > 40. & abs(eta) < 2.4'
#tight6JetCut  = 'pt > 30. & abs(eta) < 2.4'
tightJetCut   = 'pt > 30. & abs(eta) < 2.4'

## define tight JetIDs
tightCaloJetID  = '& (((jecFactor("Uncorrected") * pt) > 25. & jetID.fHPD < 0.95) | (jecFactor("Uncorrected") * pt) <= 25.) & '
tightCaloJetID += '((abs(eta) > 1.0 & (jecFactor("Uncorrected") * pt) > 55. & emEnergyFraction < 1.0) | abs(eta) <= 1.0 | (jecFactor("Uncorrected") * pt) <= 55.)'
tightPFJetID    = '& neutralHadronEnergyFraction < 0.9 & '
tightPFJetID   += 'neutralEmEnergyFraction < 0.9'

## setup the jet selection collection
tight4LeadingJets = selectedPatJets.clone(src = 'goodJetsPF',
                                          cut = tight4JetCut #+ tightPFJetID
                                          )
tight5LeadingJets = selectedPatJets.clone(src = 'goodJetsPF',
                                          cut = tight5JetCut #+ tightPFJetID
                                          )
#tightLeadingJets = selectedPatJets.clone(src = 'goodJetsPF',
#                                          cut = tight6JetCut #+ tightPFJetID
#                                          )
tight6LeadingJets = selectedPatJets.clone(src = 'goodJetsPF',
                                          cut = tight6JetCut #+ tightPFJetID
                                          )
tightLeadingJets = selectedPatJets.clone(src = 'goodJetsPF',
                                          cut = tightJetCut #+ tightPFJetID
                                          )

trackCountingHighEffBJets.src = 'goodJetsPF'
#tightBottomJets  = selectedPatJets.clone(src = 'trackCountingHighEffBJets',
#                                         cut = tight6JetCut #+ tightPFJetID
#                                         )
tightBottomJets  = selectedPatJets.clone(src = 'trackCountingHighEffBJets',
                                         cut = tightJetCut #+ tightPFJetID
                                         )

## setting up the collections for the fully-hadronic
## event selection; on these collection monitoring
## can still be performed
createJetCollections = cms.Sequence(goodJetsPF *
                                    trackCountingHighEffBJets *
                                    tight4LeadingJets *
                                    tight5LeadingJets *
                                    tight6LeadingJets *
                                    tightLeadingJets *
                                    tightBottomJets
                                    )

## ---
##    configure the cutflow scenario
## ---

## jet quality analyzer
from TopAnalysis.TopAnalyzer.JetQuality_cfi import *
## muon quality analyzer
#from TopAnalysis.TopAnalyzer.MuonQuality_cfi import *
## jet kinematics analyzer
from TopAnalysis.TopAnalyzer.JetKinematics_cfi import *
## MET kinematics analyzer
from TopAnalysis.TopAnalyzer.METKinematics_cfi import *
## muon kinematics analyzer
#from TopAnalysis.TopAnalyzer.MuonKinematics_cfi import *
## event shape analyzer
from TopAnalysis.TopAnalyzer.EventShapes_cfi import *
## genParticle analyzer
from TopAnalysis.TopAnalyzer.GenParticle_cfi import *
## ptHat analyzer
from TopAnalysis.TopAnalyzer.PtHat_cfi import *
## analyzer for special variables in fully hadronic channel
from TopAnalysis.TopAnalyzer.FullHadSpecial_cfi import *
## kinfit analyzer
#from TopAnalysis.TopAnalyzer.KinFitQuality_cfi import *
from TopAnalysis.TopAnalyzer.KinFitImprover_cfi import *
from TopAnalysis.TopAnalyzer.METKinFitAnalyzer_cfi import *
## analyzer for fully hadronic event reco
from TopAnalysis.TopAnalyzer.FullHadTopReco_cfi import *
## high level trigger filter
from TopAnalysis.TopFilter.sequences.triggerFilter_cff import *
from TopAnalysis.TopFilter.filters.NewTriggerTestFilter_cfi import *
## generator matching
from TopAnalysis.TopFilter.sequences.generatorMatching_cff import *
## kinFit producer
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff import *
from TopAnalysis.TopAnalyzer.FullHadHypothesisAnalyzer_cff import *
## PDF uncertainties
from TopAnalysis.TopAnalyzer.PDFUncertainty_cfi import *

## ---
##    TRIGGER DEFINITION
## ---

### create trigger selection for MC
#trigger = hltQuadJet30.clone( HLTPaths = ["HLT_QuadJet15U"],
#                              TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI") )
#
### create personal QuadJet40 trigger
#
## get trigger module to be added to the path
from TopAnalysis.TopUtils.patTriggerEvent_cff import *
#
### the QuadJet40 trigger itself
#hltQJ40 = filterTrigger.clone()
#
#trigger = cms.Sequence(patTriggerDefaultSequence *
#                       hltQJ40
#                       )
#
### create personal QuadJet25U trigger
#
### the QuadJet25U trigger itself
#hltQJ25U = filterTrigger.clone( whichTrigger="QuadJet25U" )
#
#trigger = cms.Sequence(patTriggerDefaultSequence *
#                       hltQJ25U
#                       )
#
### high level trigger filter
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#trigger = hltHighLevel.clone(HLTPaths = ["HLT_QuadJet25U", "HLT_QuadJet25U_v2", "HLT_QuadJet25U_v3"],
#                             throw = False)

## high level trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
trigger = hltHighLevel.clone( HLTPaths = [ #2010 trigger ('v*' to be immune to version changes)
                                           'HLT_QuadJet15U'   , 'HLT_QuadJet20U'   ,'HLT_QuadJet25U'   
                                         , 'HLT_QuadJet15U_v*', 'HLT_QuadJet20U_v*','HLT_QuadJet25U_v*'
                                         #2011 1E33 trigger ('v*' to be immune to version changes)
                                         , 'HLT_QuadJet50_Jet40_v*'
                                         #2011 1E33-2E33 trigger ('v*' to be immune to version changes)
                                         , 'HLT_QuadJet50_Jet40_Jet30_v*'
                                         #2012 5E33 trigger ('v*' to be immune to version changes)
                                         , 'HLT_QuadJet50_v*']
                                         #, 'HLT_QuadJet60_DiJet20_v*'
                                         #, 'HLT_QuadJet60_DiJet20_L1FastJet_v*']
                             , throw = False)

## ---
##    FILTER STEP 0
## ---

### vertex filter
#vertex = cms.EDFilter( "VertexSelector"
#                     , src = cms.InputTag("offlinePrimaryVertices")
#                     , cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2")
#                     , filter = cms.bool(True)
#                     )
#
### scraping filter
#noscraping = cms.EDFilter( "FilterOutScraping"
#                         , applyfilter = cms.untracked.bool(True)
#                         , debugOn = cms.untracked.bool(False)
#                         , numtrack = cms.untracked.uint32(10)
#                         , thresh = cms.untracked.double(0.25)
#                         )
#
### setup good jet selection collection
#goodJetSelection = countPatJets.clone( src = 'goodJetsPF'
#                                     , minNumber = 6
#                                     )
#
#filterStep0 = cms.Sequence( vertex
#                          * noscraping
#                          * goodJetSelection
#                          )

## to switch verbosity modes of the kinFit
#ttFullHadEvent.verbosity = 3

## configuration of kinematic fit
kinFitTtFullHadEventHypothesis.maxNComb = -1

kinFitTtFullHadEventHypothesis.bTags = 2
kinFitTtFullHadEventHypothesis.bTagAlgo = 'trackCountingHighEffBJetTags'
kinFitTtFullHadEventHypothesis.minBTagValueBJet    = 3.3
kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet = 10.2

#setForAllTtFullHadHypotheses(process, 'maxNJets', -1)
kinFitTtFullHadEventHypothesis.maxNJets = -1
ttFullHadJetPartonMatch.maxNJets        = -1

#setForAllTtFullHadHypotheses(process, 'jetCorrectionLevel', 'had')
kinFitTtFullHadEventHypothesis.jetCorrectionLevel = 'L3Absolute'
ttFullHadHypGenMatch.jetCorrectionLevel           = 'L3Absolute'
#kinFitTtFullHadEventHypothesis.jetCorrectionLevel = 'L7Parton'
#ttFullHadHypGenMatch.jetCorrectionLevel           = 'L7Parton'

#setForAllTtFullHadHypotheses(process, 'jets', 'tightLeadingJets')
kinFitTtFullHadEventHypothesis.jets = 'tightLeadingJets'
ttFullHadJetPartonMatch.jets        = 'tightLeadingJets'
ttFullHadHypGenMatch.jets           = 'tightLeadingJets'
ttFullHadHypKinFit.jets             = 'tightLeadingJets'

## exchange resolutions for jets
#from TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_cff import *
from TopAnalysis.Configuration.stringResolutions_etEtaPhi_Fall11_cff import *
kinFitTtFullHadEventHypothesis.udscResolutions = udscResolutionPF.functions
kinFitTtFullHadEventHypothesis.bResolutions    = bjetResolutionPF.functions

## configure genMatch
ttFullHadJetPartonMatch.useMaxDist = True
ttFullHadJetPartonMatch.algorithm  = 'unambiguousOnly' #'totalMinDist' #

## define ordered jets
udsall =cms.PSet(index=cms.int32(-1), correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds0   =cms.PSet(index=cms.int32(0) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds1   =cms.PSet(index=cms.int32(1) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds2   =cms.PSet(index=cms.int32(2) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds3   =cms.PSet(index=cms.int32(3) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds4   =cms.PSet(index=cms.int32(4) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
uds5   =cms.PSet(index=cms.int32(5) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("uds")   , useTree=cms.bool(False))
ball   =cms.PSet(index=cms.int32(-1), correctionLevel=cms.string('L3Absolute'), flavor=cms.string("bottom"), useTree=cms.bool(False))
bottom0=cms.PSet(index=cms.int32(0) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("bottom"), useTree=cms.bool(False))
bottom1=cms.PSet(index=cms.int32(1) , correctionLevel=cms.string('L3Absolute'), flavor=cms.string("bottom"), useTree=cms.bool(False))

## ---
##    create a monitoring sequence with name *monitoringSequence[suffix]*
## ---

listOfMonitoringSuffixes = []

def createMonitoringSequence(suffix, jetSrc, bJetSrc, state=0):

    globals()['listOfMonitoringSuffixes'] += [suffix]

    ## JET KINEMATICS

    ## kinematics analyzers
    globals()['tightBottomJetKinematics'  + suffix] = analyzeJetKinematics.clone(src = bJetSrc, analyze = ball   )
    globals()['tightLeadingJetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = udsall )
    globals()['tightLead_0_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds0   )
    globals()['tightLead_1_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds1   )
    globals()['tightLead_2_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds2   )
    globals()['tightLead_3_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds3   )
    globals()['tightLead_4_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds4   )
    globals()['tightLead_5_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = jetSrc , analyze = uds5   )
    globals()['tightBJet_0_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = bJetSrc, analyze = bottom0)
    globals()['tightBJet_1_JetKinematics' + suffix] = analyzeJetKinematics.clone(src = bJetSrc, analyze = bottom1)
    globals()['METKinematics'             + suffix] = analyzeMETKinematics.clone(srcA = 'patMETsPF')

    ## collect kinematics analyzers
    globals()['monitorKinematics' + suffix] = cms.Sequence( globals()['tightBottomJetKinematics'  + suffix]
                                                          * globals()['tightBJet_0_JetKinematics' + suffix]
                                                          * globals()['tightBJet_1_JetKinematics' + suffix]
                                                          * globals()['tightLeadingJetKinematics' + suffix]
                                                          * globals()['tightLead_0_JetKinematics' + suffix]
                                                          * globals()['tightLead_1_JetKinematics' + suffix]
                                                          * globals()['tightLead_2_JetKinematics' + suffix]
                                                          * globals()['tightLead_3_JetKinematics' + suffix]
                                                          * globals()['tightLead_4_JetKinematics' + suffix]
                                                          * globals()['tightLead_5_JetKinematics' + suffix]
                                                          * globals()['METKinematics'             + suffix]
                                                          )

    ## JET QUALITY

    ## quality analyzers
    globals()['tightBottomJetQuality'  + suffix] = analyzeJetQuality.clone(src = bJetSrc, analyze = ball   )
    globals()['tightLeadingJetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = udsall )
    globals()['tightLead_0_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds0   )
    globals()['tightLead_1_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds1   )
    globals()['tightLead_2_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds2   )
    globals()['tightLead_3_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds3   )
    globals()['tightLead_4_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds4   )
    globals()['tightLead_5_JetQuality' + suffix] = analyzeJetQuality.clone(src = jetSrc , analyze = uds5   )
    globals()['tightBJet_0_JetQuality' + suffix] = analyzeJetQuality.clone(src = bJetSrc, analyze = bottom0)
    globals()['tightBJet_1_JetQuality' + suffix] = analyzeJetQuality.clone(src = bJetSrc, analyze = bottom1)

    ## collect quality analyzers
    globals()['monitorJetsQuality' + suffix] = cms.Sequence( globals()['tightBottomJetQuality'  + suffix]
                                                           * globals()['tightBJet_0_JetQuality' + suffix]
                                                           * globals()['tightBJet_1_JetQuality' + suffix]
                                                           * globals()['tightLeadingJetQuality' + suffix]
                                                           * globals()['tightLead_0_JetQuality' + suffix]
                                                           * globals()['tightLead_1_JetQuality' + suffix]
                                                           * globals()['tightLead_2_JetQuality' + suffix]
                                                           * globals()['tightLead_3_JetQuality' + suffix]
                                                           * globals()['tightLead_4_JetQuality' + suffix]
                                                           * globals()['tightLead_5_JetQuality' + suffix]
                                                           )

    ## EVENT SHAPES

    ## collect event shape analyzers
    globals()['eventShapes' + suffix] = analyzeEventShapes.clone(src = jetSrc)

    ## monitor sequence for event shape analyzers
    globals()['monitorEventShapes' + suffix] = cms.Sequence( globals()['eventShapes' + suffix] )

    ## FULL HAD SPECIAL

    ## collect analyzers specially for full hadronic analysis
    globals()['fullHadSpecial' + suffix] = analyzeFullHadSpecials.clone(src = jetSrc)

    ## monitor sequence for specially for full hadronic analyzers
    globals()['monitorFullHadSpecials' + suffix] = cms.Sequence( globals()['fullHadSpecial' + suffix] )

    ## GENERATOR INFORMATION

    ## collect analyzers for genParticles
    globals()['genParticles' + suffix] = analyzeGenParticles.clone()

    ## collect analyzers for ptHat
    globals()['ptHat' + suffix] = analyzePtHat.clone()

    ## monitor sequence for generator
    globals()['monitorGenerator' + suffix] = cms.Sequence( globals()['genParticles' + suffix]
                                                         * globals()['ptHat' + suffix]
                                                         )

    ## collect all monitoring sequences
    globals()['monitoringSequence' + suffix] = cms.Sequence( globals()['monitorKinematics'      + suffix]
                                                           * globals()['monitorJetsQuality'     + suffix]
                                                           * globals()['monitorEventShapes'     + suffix]
                                                           * globals()['monitorFullHadSpecials' + suffix]
                                                           * globals()['monitorGenerator'       + suffix]
                                                           )

    if(state > 0):
        
        ## KINFIT QUALITY AND FULL HAD TOP RECO

        ## kinfit quality analyzer
        ## collect kinfit quality analyzers
        #globals()['kinFitQuality'   + suffix] = analyzeKinFitQuality.clone (srcB = jetSrc)
        globals()['kinFitImprover0' + suffix] = analyzeKinFitImprover.clone(srcB = jetSrc, analyze = cms.PSet(comboType = cms.uint32(0) ) )
        globals()['kinFitImprover1' + suffix] = analyzeKinFitImprover.clone(srcB = jetSrc, analyze = cms.PSet(comboType = cms.uint32(1) ) )
        globals()['kinFitImprover2' + suffix] = analyzeKinFitImprover.clone(srcB = jetSrc, analyze = cms.PSet(comboType = cms.uint32(2) ) )
        globals()['kinFitImprover3' + suffix] = analyzeKinFitImprover.clone(srcB = jetSrc, analyze = cms.PSet(comboType = cms.uint32(3) ) )
        globals()['kinFitImprover4' + suffix] = analyzeKinFitImprover.clone(srcB = jetSrc, analyze = cms.PSet(comboType = cms.uint32(4) ) )

        #globals()['kinFitQuality' + suffix].analyze.udscResolutions = udscResolutionPF.functions
        #globals()['kinFitQuality' + suffix].analyze.bResolutions    = bjetResolutionPF.functions
        #globals()['kinFitQuality' + suffix].analyze.udscResolutions = udscResolution.functions
        #globals()['kinFitQuality' + suffix].analyze.bResolutions    = bjetResolution.functions

        ## collect fully hadronic top reco analyzers
        globals()['fullHadTopReco' + suffix] = analyzeFullHadTopReco.clone(srcB = jetSrc)
        globals()['fullHadTopReco' + suffix].analyze.bTagAlgo = 'trackCountingHighEffBJetTags' 
        globals()['METKinFit'      + suffix] = analyzeMETKinFit.clone(JetSrc = jetSrc)
        
        ## monitor sequence for kinfit quality analyzers
        globals()['monitorKinFit' + suffix] = cms.Sequence(# globals()['kinFitQuality'   + suffix]
                                                            globals()['kinFitImprover0' + suffix]
                                                          * globals()['kinFitImprover1' + suffix]
                                                          * globals()['kinFitImprover2' + suffix]
                                                          * globals()['kinFitImprover3' + suffix]
                                                          * globals()['kinFitImprover4' + suffix]
                                                          * globals()['METKinFit'       + suffix]
                                                          * globals()['fullHadTopReco'  + suffix]
                                                          )

        globals()['monitoringSequence' + suffix] = cms.Sequence( globals()['monitorKinematics'      + suffix]
                                                               * globals()['monitorJetsQuality'     + suffix]
                                                               * globals()['monitorEventShapes'     + suffix]
                                                               * globals()['monitorFullHadSpecials' + suffix]
                                                               * globals()['monitorGenerator'       + suffix]
                                                               * globals()['monitorKinFit' + suffix]
                                                               )

    ## add monitoring sequence to full analysis sequence
    globals()['analyseFullHadronicSelection'] += globals()['monitoringSequence' + suffix]



## ---
##    MONITOR STEP 0
## ---

## PDF WEIGHTS for PDF uncertainties
PDFUncertainty_0 = analyzePDFUncertainty.clone()

## ---
##    Final Sequence to be Run
##
##    further elements added *in time*, see below
## ---
analyseFullHadronicSelection = cms.Sequence( PDFUncertainty_0
                                           ## do the hlt triggering
                                           * trigger
                                           ## create the jet collections
                                           * createJetCollections
                                           ## do the matching
                                           #* matchJetsToPartons
                                           ## filter step 0
                                           #* filterStep0
                                           )

createMonitoringSequence("_0", 'goodJetsPF', 'trackCountingHighEffBJets')

## ---
##    FILTER STEP 1
##    Multi-Jet Selection
## ---

## select all events with at least 4 jets
leading4JetSelection = countPatJets.clone( src = 'tight4LeadingJets'
                                         , minNumber = 4
                                         )

## select all events with at least 5 jets
leading5JetSelection = countPatJets.clone( src = 'tight5LeadingJets'
                                         , minNumber = 5
                                         )

## select all events with at least 6 jets
#leading6JetSelection = countPatJets.clone( src = 'tightLeadingJets'
#                                         , minNumber = 6
#                                         )
leading6JetSelection = countPatJets.clone( src = 'tight6LeadingJets'
                                         , minNumber = 6
                                         )

## select all events with at least 6 jets
leadingJetSelection = cms.Sequence( leading4JetSelection
                                  * leading5JetSelection
                                  * leading6JetSelection
                                  )

analyseFullHadronicSelection += leadingJetSelection

## ---
##    MONITOR STEP 1
##    Multi-Jet Selection
## ---

createMonitoringSequence("_1", 'tightLeadingJets', 'tightBottomJets')

## ---
##    FILTER STEP 2
##    Multi-Jet MVA Selection
## ---

from TopQuarkAnalysis.TopEventSelection.TtFullHadSignalSelMVAComputer_cff import *
from TopAnalysis.TopFilter.filters.MVADiscriminatorFilter_cfi import *

mvaDisc = findTtFullHadSignalSelMVA.clone(jets = 'tightLeadingJets')
filterMVADiscriminator.maxDisc = 0.92

TtFullHadSignalSelMVAFileSource.ttFullHadSignalSelMVA = 'TopAnalysis/Configuration/data/TtFullHadSignalSel_tmva_test.mva'
mvaSelection = cms.Sequence( mvaDisc
                           * filterMVADiscriminator
                           )

analyseFullHadronicSelection += mvaSelection

#dump = cms.EDAnalyzer('EventContentAnalyzer')
#analyseFullHadronicSelection += dump

## ---
##    MONITOR STEP 2
##    Multi-Jet MVA Selection
## ---

createMonitoringSequence("_2", 'tightLeadingJets', 'tightBottomJets')

## ---
##    FILTER STEP 3
##    B-Tag Selection
## ---

## select events with at least 2 b jets
tightBottomJetSelection   = countPatJets.clone( src = 'tightBottomJets'
                                              , minNumber = 2
                                              )

bottomJetSelection = cms.Sequence(tightBottomJetSelection
                                  )

analyseFullHadronicSelection += bottomJetSelection

## do kinematic fit on all events with 2 b-tags
analyseFullHadronicSelection += makeTtFullHadEvent

## ---
##    MONITOR STEP 3
##    B-Tag Selection
## ---

createMonitoringSequence("_3", 'tightLeadingJets', 'tightBottomJets', 1)

## ---
##    FILTER STEP 4
##    Kinematic Fit
## ---

## kinfit quality filter
from TopQuarkAnalysis.TopEventProducers.producers.TtFullHadEvtFilter_cfi import *
filterKinFitQuality = ttFullHadEventFilter.clone( cut = cms.string("isHypoValid('kKinFit') & fitProb > 0.01") )
#filterKinFitQuality = ttFullHadEventFilter.clone( cut = cms.string("isHypoValid('kKinFit') & fitChi2 < 30.0") )

analyseFullHadronicSelection += filterKinFitQuality

## ---
##    MONITOR STEP 4
##    Kinematic Fit
## ---

createMonitoringSequence("_4", 'tightLeadingJets', 'tightBottomJets', 1)

## ---
##    FILTER STEP 5 (not used)
## ---

## event shape filter
from TopAnalysis.TopFilter.filters.EventShapeFilter_cfi import *
filterEventShapes = filterEventShape.clone( minC = 0.75 )

## ---
##    MONITOR STEP 5
## ---

## To be added in time



## PDF WEIGHTS for PDF uncertainties
PDFUncertainty_3 = analyzePDFUncertainty.clone()
analyseFullHadronicSelection += PDFUncertainty_3

## ---
##    provide a function to use trees if available instead of histograms
## ---
def useTreesAsOutput(process):
    from PhysicsTools.PatAlgos.tools.helpers import listModules
    for mod in listModules(process.analyseFullHadronicSelection):
        if(hasattr(mod, 'analyze')):
            if(hasattr(mod.analyze, 'useTree')):
                mod.analyze.useTree = True
            if(hasattr(mod, 'useTree')):
                mod.useTree = True

## ---
##    provide a function to disable parts of the selection
## ---
def disableCountFilter(whichCountFilter):
    whichCountFilter.minNumber = 0
    whichCountFilter.maxNumber = 999999

## ---
##    run on data
## ---
def runOnData(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing all elements from the sequence '
    print '*analyseFullHadronicSelection* that rely '
    print 'on generator information to run properly '
    print '++++++++++++++++++++++++++++++++++++++++++++'

    ## removal of generator based information
    process.analyseFullHadronicSelection.remove(process.matchJetsToPartons)

    for suf in listOfMonitoringSuffixes:
        process.analyseFullHadronicSelection.remove(getattr(process,'monitorGenerator'+suf))

    ## switch to residual jet energy correction for data
    #if(hasattr(process, 'goodJetsPF')):
    process.kinFitTtFullHadEventHypothesis.jetCorrectionLevel = 'L2L3Residual'
    process.ttFullHadHypGenMatch.jetCorrectionLevel           = 'L2L3Residual'
    if(hasattr(process, 'analyzeFullHadEventMixer')):
        process.analyzeFullHadEventMixer.jetCorrectionLevel   = 'L2L3Residual'
    
    udsall.correctionLevel  = 'L2L3Residual'
    uds0.correctionLevel    = 'L2L3Residual'
    uds1.correctionLevel    = 'L2L3Residual'
    uds2.correctionLevel    = 'L2L3Residual'
    uds3.correctionLevel    = 'L2L3Residual'
    uds4.correctionLevel    = 'L2L3Residual'
    uds5.correctionLevel    = 'L2L3Residual'
    ball.correctionLevel    = 'L2L3Residual'
    bottom0.correctionLevel = 'L2L3Residual'
    bottom1.correctionLevel = 'L2L3Residual'

    ## changes needed for analysis
    #trigger.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
    #process.patTrigger.processName = 'HLT'
    #vertex.cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2")

    ## to switch verbosity modes of the kinFit
    #process.ttFullHadEvent.verbosity = 3


## ---
##    remove modules that produce monitoring plots during the cutflow
## ---
def removeMonitoringOfCutflow(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing all monitoring elements from the '
    print 'sequence *analyseFullHadronicSelection* so'
    print 'only a pure selection of events is done '
    print '++++++++++++++++++++++++++++++++++++++++++++'
    for suf in listOfMonitoringSuffixes:
        process.analyseFullHadronicSelection.remove(getattr(process,'monitoringSequence'+suf))

    process.analyseFullHadronicSelection.remove(process.PDFUncertainty_0)
    process.analyseFullHadronicSelection.remove(process.PDFUncertainty_3)

## ---
##    remove default trigger
## ---
def removeTrigger(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing the trigger from standard '
    print 'fully hadronic event selection '
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.remove(process.trigger)
    
## ---
##    switch all necessary filters to run this sequence for background estimation
## ---
def runAsBackgroundEstimation(process, whichEstimate):

    if whichEstimate == 0 :
        print '++++++++++++++++++++++++++++++++++++++++++++'
        print 'no background estimate done'
        print '++++++++++++++++++++++++++++++++++++++++++++'
    elif whichEstimate == 1 :
        print '++++++++++++++++++++++++++++++++++++++++++++'
        print 'loading background estimation and switching'
        print 'bottomJetSelection, implemation into'
        print 'sequence needs to be done BY HAND !!!'
        print '++++++++++++++++++++++++++++++++++++++++++++'
        process.tightBottomJetSelection.minNumber = 0
        process.tightBottomJetSelection.maxNumber = 0
    
        ## ---
        ##    load FullHadQCDEstimation Analyzer
        ## ---
        
        ## analyze the QCD estimation for fully hadronic ttbar events
        process.load("TopAnalysis.TopAnalyzer.FullHadQCDEstimation_cfi")
        process.analyzeFullHadQCDEstimation = process.analyzeFullHadQCDEstimation.clone(JetSrc = "tightLeadingJets")
        #process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHPT40"
        #process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHEM40"
        #process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHEM40MC"
        process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHEM30"
        #process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHEM30clean"
        
        #process.load("TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_cff")
        process.analyzeFullHadQCDEstimation.udscResolutions = udscResolutionPF.functions
        process.analyzeFullHadQCDEstimation.bResolutions    = bjetResolutionPF.functions
        #process.analyzeFullHadQCDEstimation.udscResolutions = udscResolution.functions
        #process.analyzeFullHadQCDEstimation.bResolutions    = bjetResolution.functions

        process.analyseFullHadronicSelection += process.analyzeFullHadQCDEstimation

    elif whichEstimate == 2 :
        ## ---
        ##    load FullHadEventMixer
        ## ---
       
        print '++++++++++++++++++++++++++++++++++++++++++++'
        print 'loading background estimation using the '
        print 'event mixing tequnique, implemation into'
        print 'sequence needs to be done BY HAND !!!'
        print '++++++++++++++++++++++++++++++++++++++++++++'
        
        process.load("TopAnalysis.TopAnalyzer.FullHadEventMixer_cfi")
        #process.analyzeFullHadEventMixer = process.analyzeFullHadEventMixer.clone(JetSrc = "tightLeadingJets", jetCorrectionLevel = 'L7Parton')
        process.analyzeFullHadEventMixer = process.analyzeFullHadEventMixer.clone(JetSrc = "tightLeadingJets")
        
        #from TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_cff import *
        process.analyzeFullHadEventMixer.udscResolutions = udscResolutionPF.functions
        process.analyzeFullHadEventMixer.bResolutions    = bjetResolutionPF.functions

        process.analyseFullHadronicSelection += process.analyzeFullHadEventMixer

    else:
        print 'whichEstimate =', whichEstimate, 'not allowed, only supported options: 0 (no background estimate), 1 (b-tag weighting), 2 (event mixing)'


## ---
##    run analysis on CaloJets & MET instead of PFJets & PFMet
## ---
def runOnCalo(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'switching all inputs in to run on Calo'
    print 'goodJetsPF -> goodJets'
    print 'patMETsPF  -> patMETs'
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.replace(process.goodJetsPF, process.goodJets)

    ## exchange PFJetID cuts to JetID
    process.tight4LeadingJets.cut = tight4JetCut #+ tightCaloJetID
    process.tight5LeadingJets.cut = tight5JetCut #+ tightCaloJetID
    #process.tightLeadingJets.cut  = tight6JetCut #+ tightCaloJetID
    process.tight6LeadingJets.cut = tight6JetCut #+ tightCaloJetID
    process.tightLeadingJets.cut  = tightJetCut #+ tightCaloJetID
    #process.tightBottomJets.cut   = tight6JetCut #+ tightCaloJetID
    process.tightBottomJets.cut   = tightJetCut  #+ tightCaloJetID

    ## exchange resolutions for CaloJets
    from TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_cff import udscResolution,bjetResolution
    process.kinFitTtFullHadEventHypothesis.udscResolutions = udscResolution.functions
    process.kinFitTtFullHadEventHypothesis.bResolutions    = bjetResolution.functions

    #for suf in listOfMonitoringSuffixes:
    #    if(hasattr(process, 'kinFitQuality'+suf)):
    #        getattr(process, 'kinFitQuality'+suf).analyze.udscResolutions = udscResolution.functions
    #        getattr(process, 'kinFitQuality'+suf).analyze.bResolutions    = bjetResolution.functions

    ## run kinematic fit for CaloJets with L1L2L3(Res)L5 correted jets
    process.kinFitTtFullHadEventHypothesis.jetCorrectionLevel = 'L5Hadron'
    process.ttFullHadHypGenMatch.jetCorrectionLevel           = 'L5Hadron'

    ## replace PFJets and PFMET with CaloJets and CaloMET
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.analyseFullHadronicSelection, 'goodJetsPF', 'goodJets')
    massSearchReplaceAnyInputTag(process.analyseFullHadronicSelection, 'patMETsPF' , 'patMETs')

    ## give correct resolutions for background estimation if done
    if(hasattr(process, 'analyzeFullHadQCDEstimation')):
        process.analyzeFullHadQCDEstimation.analyze.udscResolutions = udscResolution.functions
        process.analyzeFullHadQCDEstimation.analyze.bResolutions    = bjetResolution.functions

    ## in case of mc, smear the energy resolution additionally
    if(hasattr(process, 'scaledJetEnergy')):
        process.scaledJetEnergy.inputJets = "selectedPatJets"
        process.scaledJetEnergy.inputMETs = "patMETs"
        process.scaledJetEnergy.payload   = "AK5Calo"

## ---
##    modify bTagging discriminator for bTaggers
## ---
def modifyBTagDiscs(process, algo, newMinDisc, newMaxDisc):

    getattr(process, algo + "BJets").cut = 'bDiscriminator("' + algo + 'BJetTags") > ' + str(newMinDisc)
    process.analyseFullHadronicSelection.replace(process.trackCountingHighEffBJets, getattr(process, algo + "BJets"))

    for suf in listOfMonitoringSuffixes:
        if(hasattr(process, 'fullHadTopReco'+suf)):
            getattr(process, 'fullHadTopReco'+suf).analyze.bTagAlgo = algo + 'BJetTags' 

    process.kinFitTtFullHadEventHypothesis.bTagAlgo            = algo + 'BJetTags'
    process.kinFitTtFullHadEventHypothesis.minBTagValueBJet    = newMinDisc
    process.kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet = newMaxDisc

    if(hasattr(process, 'analyzeFullHadEventMixer')):
        process.analyzeFullHadEventMixer.bTagAlgo            = algo + 'BJetTags' 
        process.analyzeFullHadEventMixer.minBTagValueBJet    = newMinDisc
        process.analyzeFullHadEventMixer.maxBTagValueNonBJet = newMaxDisc

    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.analyseFullHadronicSelection, 'trackCountingHighEffBJets', algo + 'BJets')

def modifyBTagDisc(process, algo, newDisc):

    modifyBTagDiscs(process, algo, newDisc, process.kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet.value() * newDisc / process.kinFitTtFullHadEventHypothesis.minBTagValueBJet.value())

## ---
##    switch to trackCountingHighEfficiency bTagger
## ---
def switchToTCHEL(process):
    modifyBTagDiscs(process, 'trackCountingHighEff', 1.9, 3.3)

## ---
##    switch to trackCountingHighEfficiency bTagger
## ---
def switchToTCHEM(process):
    modifyBTagDiscs(process, 'trackCountingHighEff', 3.3, 10.2)

## ---
##    switch to trackCountingHighEfficiency bTagger
## ---
def switchToTCHPM(process):
    modifyBTagDiscs(process, 'trackCountingHighPur', 1.93, 3.41)

## ---
##    switch to trackCountingHighEfficiency bTagger
## ---
def switchToTCHPT(process):
    modifyBTagDiscs(process, 'trackCountingHighPur', 3.41, 3.41)

## ---
##    switch to simpleSecondaryVertexHighEff bTagger
## ---
def switchToSSVHEM(process):
    modifyBTagDiscs(process, 'simpleSecondaryVertexHighEff', 1.74, 3.05)

## ---
##    switch to simpleSecondaryVertexHighPur bTagger
## ---
def switchToSSVHPT(process):
    modifyBTagDiscs(process, 'simpleSecondaryVertexHighPur', 2.00, 2.00)

## ---
##    switch to combinedSecondaryVertex bTagger
## ---
def switchToCSVL(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertex', 0.244, 0.679)

## ---
##    switch to combinedSecondaryVertex bTagger
## ---
def switchToCSVM(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertex', 0.679, 0.898)

## ---
##    switch to combinedSecondaryVertex bTagger
## ---
def switchToCSVT(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertex', 0.898, 0.898)
def switchToCSVTSysUp(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertex', 0.891, 0.891)
def switchToCSVTSysDown(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertex', 0.904, 0.904)

## ---
##    switch to combinedSecondaryVertexMVA bTagger
## ---
def switchToCSVMVA(process):
    modifyBTagDiscs(process, 'combinedSecondaryVertexMVA', 0.75, 0.921) # self-derived WORKING POINTS

## ---
##    increase resolutions of kinematic fit
## ---
def increaseKinFitResolution(process, binning, factors):
    #process.kinFitTtFullHadEventHypothesis.energyResolutionSmearFactor = factor
    process.kinFitTtFullHadEventHypothesis.jetEnergyResolutionScaleFactors = cms.vdouble(factors)
    process.kinFitTtFullHadEventHypothesis.jetEnergyResolutionEtaBinning   = cms.vdouble(binning)

    if(hasattr(process, 'analyzeFullHadEventMixer')):
        process.analyzeFullHadEventMixer.jetEnergyResolutionScaleFactors = cms.vdouble(factors)
        process.analyzeFullHadEventMixer.jetEnergyResolutionEtaBinning   = cms.vdouble(binning)
    if(hasattr(process, 'analyzeFullHadQCDEstimation')):
        process.analyzeFullHadQCDEstimation.jetEnergyResolutionScaleFactors = cms.vdouble(factors)
        process.analyzeFullHadQCDEstimation.jetEnergyResolutionEtaBinning   = cms.vdouble(binning)


## ---
##    remove PDF uncertainty histograms / trees
## ---
def removePDFUncertainties(process):
    process.analyseFullHadronicSelection.remove(process.PDFUncertainty_0)
    process.analyseFullHadronicSelection.remove(process.PDFUncertainty_3)

## ---
##    switch to b' search mode (pp -> b'b' -> (bZ) (bZ)
## ---
def bPrimeSearchMode(process):
    process.kinFitTtFullHadEventHypothesis.mW = 91.2
    process.kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet = process.kinFitTtFullHadEventHypothesis.minBTagValueBJet

