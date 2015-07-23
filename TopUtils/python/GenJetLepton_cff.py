import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets


genParticlesForJets = cms.EDFilter("CandViewSelector",
                       src = cms.InputTag("genParticles"),
                       cut = cms.string("!((abs(pdgId) = 11 | abs(pdgId) = 13) & pt > 20 & abs(eta) < 2.5)")
                       )

myAk5GenJetsNoLp = ak5GenJets.clone( src = cms.InputTag("genParticlesForJets") )


recoGenJets = cms.Sequence(ak5GenJets*
                           myAk5GenJetsNoLp
                          )
