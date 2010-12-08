import FWCore.ParameterSet.Config as cms

## jet selector
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *

## ---
##    setup jet quality studies
## ---

noOverlapJetsPF = cleanPatJets.clone(
    src = cms.InputTag("selectedPatJetsAK5PF"), 

    # preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),

    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("vertexSelectedMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string('pt > 20. & abs(eta) < 2.1 &'
                                            'combinedMuon.isNull = 0 &'
                                            'isTrackerMuon() =1 &'
                                            '(trackIso+caloIso)/pt < 0.05 &'
                                            'innerTrack.numberOfValidHits >= 11 &'
                                            'globalTrack.normalizedChi2 < 10.0 &'
                                            'globalTrack.hitPattern.numberOfValidMuonHits>0 &'
                                            'abs(dB)<0.02 &'
                                            'innerTrack.hitPattern.pixelLayersWithMeasurement>=1 &'
                                            'numberOfMatches>1'),
           deltaR              = cms.double(0.1),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        )
    ),
    # finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
)

## getting started
centralJets  = selectedPatJets.clone(src = 'selectedPatJets',
                                     cut = 'abs(eta) < 2.4'
                                     )
reliableJets = selectedPatJets.clone(src = 'selectedPatJets',
                                     cut = 'abs(eta) < 2.4 & pt > 30.'
                                     )
goodJets     = selectedPatJets.clone(src = 'selectedPatJets',
                                     cut = 'abs(eta) < 2.4 & pt > 30. &'
                                           'emEnergyFraction > 0.01   &'
                                           'jetID.fHPD < 0.98         &'
                                           'jetID.n90Hits > 1'
                                     )
goodJetsPF = selectedPatJets.clone(src = 'selectedPatJetsAK5PF',
                                   cut = 'abs(eta) < 2.4 & pt > 20. &'
                                         'chargedHadronEnergyFraction > 0.0  &'
                                         'neutralHadronEnergyFraction < 0.99 &'
                                         'chargedEmEnergyFraction     < 0.99 &'
                                         'neutralEmEnergyFraction     < 0.99 &'
                                         'chargedMultiplicity > 0            &'
                                         'nConstituents > 1'
                                   )
centralJetsPF20  = selectedPatJets.clone(src = 'noOverlapJetsPF',
                                     cut = 'abs(eta) < 2.4'
                                     )
centralJetsPF25  = selectedPatJets.clone(src = 'centralJetsPF20',
                                     cut = 'pt > 25.'
                                     )
centralJetsPF30  = selectedPatJets.clone(src = 'centralJetsPF20',
                                     cut = 'pt > 30.'
                                     )
reliableJetsPF20 = selectedPatJets.clone(src = 'noOverlapJetsPF',
                                     cut = 'abs(eta) < 2.4 & pt > 20.'
                                     )
reliableJetsPF25 = selectedPatJets.clone(src = 'reliableJetsPF20',
                                     cut = 'pt > 25.'
                                     )
reliableJetsPF30 = selectedPatJets.clone(src = 'reliableJetsPF20',
                                     cut = 'pt > 30.'
                                     )
goodJetsPF20     = selectedPatJets.clone(src = 'noOverlapJetsPF',
                                         cut = 'abs(eta) < 2.4 & pt > 20.          &'
                                               'chargedHadronEnergyFraction > 0.0  &'
                                               'neutralHadronEnergyFraction < 0.99 &'
                                               'chargedEmEnergyFraction     < 0.99 &'
                                               'neutralEmEnergyFraction     < 0.99 &'
                                               'chargedMultiplicity > 0            &'
                                               'nConstituents > 1'
                                         )

goodJetsPF25     = selectedPatJets.clone(src = 'goodJetsPF20',
                                         cut = 'pt > 25.'
                                         )
goodJetsPF30     = selectedPatJets.clone(src = 'goodJetsPF20',
                                         cut = 'pt > 30.'
                                         )
## N-1 collections
noEtaJets      = selectedPatJets.clone(src = 'selectedPatJets',
                                       cut = 'pt > 30.                  &'
                                             '0.01 < emEnergyFraction   &'
                                             'jetID.fHPD < 0.98         &'
                                             'jetID.n90Hits > 1'
                                       )
noPtJets       = selectedPatJets.clone(src = 'selectedPatJets',
                                       cut = 'abs(eta) < 2.4            &'
                                             '0.01 < emEnergyFraction   &'
                                             'jetID.fHPD < 0.98         &'
                                             'jetID.n90Hits > 1'
                                       )
noEmJets       = selectedPatJets.clone(src = 'selectedPatJets',
                                       cut = 'pt > 30. & abs(eta) < 2.4 &'
                                             'jetID.fHPD < 0.98         &'
                                             'jetID.n90Hits > 1'
                                       )
noN90HitsJets  = selectedPatJets.clone(src = 'selectedPatJets',
                                       cut = 'pt > 30. & abs(eta) < 2.4 &'
                                             '0.01 < emEnergyFraction   &'
                                             'jetID.fHPD < 0.98          '                                      
                                       )
nofHPDJets     = selectedPatJets.clone(src = 'selectedPatJets',
                                       cut = 'pt > 30. & abs(eta) < 2.4 &'
                                             '0.01 < emEnergyFraction   &'
                                             'jetID.n90Hits > 1'                                 
                                       )

## a kinematically well defined jet with
## reliable calibration and a robust rej
## of photons, electrons and prompt pi0's
selectGoodJets = cms.Sequence(noOverlapJetsPF  *
                              centralJets      *
                              reliableJets     *
                              goodJets         *
                              goodJetsPF       *
                              centralJetsPF20  *
                              reliableJetsPF20 *
                              goodJetsPF20     *
                              centralJetsPF25  *
                              reliableJetsPF25 *
                              goodJetsPF25     *
                              centralJetsPF30  *
                              reliableJetsPF30 *
                              goodJetsPF30
                              )

## collect the N-1 collections
selectNMinusOneJets = cms.Sequence(noEtaJets*
                                   noPtJets*
                                   noEmJets*
                                   noN90HitsJets*
                                   nofHPDJets
                                   )


## check for different btag properties
trackCountingHighPurBJets       = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"trackCountingHighPurBJetTags\") > 1.93'
                                                        )
trackCountingHighEffBJets       = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"trackCountingHighEffBJetTags\") > 3.3'
                                                        )
simpleSecondaryVertexBJets      = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"simpleSecondaryVertexBJetTags\") > 1.74'
                                                        )
simpleSecondaryVertexNegBJets   = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"simpleSecondaryVertexNegativeBJetTags\") > 3.0'
                                                        )
combinedSecondaryVertexBJets    = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"combinedSecondaryVertexBJetTags\") > 0.8'
                                                        )
combinedSecondaryVertexMVABJets = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") > 0.575'
                                                        )
softMuonBJets                   = selectedPatJets.clone(src = 'goodJets',
                                                        cut = 'bDiscriminator(\"softMuonBJetTags\") > 0.3'
                                                        )

## a goodJet fullfilling different btag
## criteria
selectBTaggedJets = cms.Sequence(goodJets                      *
                                 trackCountingHighPurBJets     *
                                 trackCountingHighEffBJets     *
                                 simpleSecondaryVertexBJets    *
                                 simpleSecondaryVertexNegBJets *
                                 combinedSecondaryVertexBJets  *
                                 combinedSecondaryVertexMVABJets  *
                                 softMuonBJets
                                 )
