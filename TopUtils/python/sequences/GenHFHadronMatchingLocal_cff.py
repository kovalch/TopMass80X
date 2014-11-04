import FWCore.ParameterSet.Config as cms


from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from TopAnalysis.TopUtils.GenJetParticles_cfi import genParticlesForJetsPlusNoHadron
from TopAnalysis.TopUtils.GenJetParticles_cff import genParticlesForJetsNoNuPlusNoHadron
from TopAnalysis.TopUtils.GenHFHadronMatcherLocal_cfi import matchGenHFHadronLocal


# Supplies PDG ID to real name resolution of MC particles
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *




# Configuration for matching B-hadrons ================================================================
genParticlesForJetsPlusBHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusBHadron.injectHadronFlavours = [5]
ak5GenJetsPlusBHadron = ak5GenJets.clone()
ak5GenJetsPlusBHadron.src = "genParticlesForJetsPlusBHadron"
matchGenBHadronLocal = matchGenHFHadronLocal.clone()
matchGenBHadronLocal.flavour = 5
matchGenBHadronLocal.genJets = cms.InputTag("ak5GenJetsPlusBHadron", "", "")
genBHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusBHadron * ak5GenJetsPlusBHadron * matchGenBHadronLocal )


# Configuration for matching C-hadrons =================================================================
genParticlesForJetsPlusCHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusCHadron.injectHadronFlavours = [4]
ak5GenJetsPlusCHadron = ak5GenJets.clone()
ak5GenJetsPlusCHadron.src = "genParticlesForJetsPlusCHadron"
matchGenCHadronLocal = matchGenHFHadronLocal.clone()
matchGenCHadronLocal.flavour = 4
matchGenCHadronLocal.genJets = cms.InputTag("ak5GenJetsPlusCHadron", "", "")
genCHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusCHadron * ak5GenJetsPlusCHadron * matchGenCHadronLocal )


# Configuration for matching B- and C-hadrons ==========================================================
genParticlesForJetsPlusBCHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusBCHadron.injectHadronFlavours = [5, 4]
ak5GenJetsPlusBCHadron = ak5GenJets.clone()
ak5GenJetsPlusBCHadron.src = "genParticlesForJetsPlusBCHadron"
matchGenBCHadronBLocal = matchGenBHadronLocal.clone()
matchGenBCHadronBLocal.genJets = cms.InputTag("ak5GenJetsPlusBCHadron", "", "")
matchGenBCHadronCLocal = matchGenCHadronLocal.clone()
matchGenBCHadronCLocal.genJets = cms.InputTag("ak5GenJetsPlusBCHadron", "", "")
genBCHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusBCHadron * ak5GenJetsPlusBCHadron * matchGenBCHadronBLocal * matchGenBCHadronCLocal )


# Configuration for matching B-hadrons (jets without neutrinos) ========================================
genParticlesForJetsNoNuPlusBHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusBHadron.injectHadronFlavours = [5]
ak5GenJetsNoNuPlusBHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusBHadron.src = "genParticlesForJetsNoNuPlusBHadron"
matchGenBHadronNoNuLocal = matchGenBHadronLocal.clone()
matchGenBHadronNoNuLocal.genJets = cms.InputTag("ak5GenJetsNoNuPlusBHadron", "", "")
genBHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusBHadron * ak5GenJetsNoNuPlusBHadron * matchGenBHadronNoNuLocal )


# Configuration for matching C-hadrons (jets without neutrinos) ========================================
genParticlesForJetsNoNuPlusCHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusCHadron.injectHadronFlavours = [4]
ak5GenJetsNoNuPlusCHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusCHadron.src = "genParticlesForJetsNoNuPlusCHadron"
matchGenCHadronNoNuLocal = matchGenCHadronLocal.clone()
matchGenCHadronNoNuLocal.genJets = cms.InputTag("ak5GenJetsNoNuPlusCHadron", "", "")
genCHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusCHadron * ak5GenJetsNoNuPlusCHadron * matchGenCHadronNoNuLocal )


# Configuration for matching B- and C-hadrons (jets without neutrinos) =================================
genParticlesForJetsNoNuPlusBCHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusBCHadron.injectHadronFlavours = [5, 4]
ak5GenJetsNoNuPlusBCHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusBCHadron.src = "genParticlesForJetsNoNuPlusBCHadron"
matchGenBCHadronBNoNuLocal = matchGenBHadronLocal.clone()
matchGenBCHadronBNoNuLocal.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
matchGenBCHadronCNoNuLocal = matchGenCHadronLocal.clone()
matchGenBCHadronCNoNuLocal.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
genBCHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusBCHadron * ak5GenJetsNoNuPlusBCHadron * matchGenBCHadronBNoNuLocal * matchGenBCHadronCNoNuLocal )


