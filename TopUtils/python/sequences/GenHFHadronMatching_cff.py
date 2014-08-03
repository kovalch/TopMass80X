import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff import ak5GenJets
from TopAnalysis.TopUtils.GenJetParticles_cfi import genParticlesForJetsPlusNoHadron
from TopAnalysis.TopUtils.GenJetParticles_cff import genParticlesForJetsNoNuPlusNoHadron
from TopAnalysis.TopUtils.GenHFHadronMatcher_cfi import matchGenHFHadron


# Configuration for matching B-hadrons ================================================================
genParticlesForJetsPlusBHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusBHadron.injectHadronFlavours = cms.vint32(5)
ak5GenJetsPlusBHadron = ak5GenJets.clone()
ak5GenJetsPlusBHadron.src = "genParticlesForJetsPlusBHadron"
matchGenBHadron = matchGenHFHadron.clone()
matchGenBHadron.flavour = 5
matchGenBHadron.genJets = cms.InputTag("ak5GenJetsPlusBHadron", "", "")
GenBHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusBHadron * ak5GenJetsPlusBHadron * matchGenBHadron )

# Configuration for matching C-hadrons =================================================================
genParticlesForJetsPlusCHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusCHadron.injectHadronFlavours = cms.vint32(4)
ak5GenJetsPlusCHadron = ak5GenJets.clone()
ak5GenJetsPlusCHadron.src = "genParticlesForJetsPlusCHadron"
matchGenCHadron = matchGenHFHadron.clone()
matchGenCHadron.flavour = 4
matchGenCHadron.genJets = cms.InputTag("ak5GenJetsPlusCHadron", "", "")
GenCHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusCHadron * ak5GenJetsPlusCHadron * matchGenCHadron )

# Configuration for matching B- and C-hadrons ==========================================================
genParticlesForJetsPlusBCHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsPlusBCHadron.injectHadronFlavours = cms.vint32(5, 4)
ak5GenJetsPlusBCHadron = ak5GenJets.clone()
ak5GenJetsPlusBCHadron.src = "genParticlesForJetsPlusBCHadron"
matchGenBCHadronB = matchGenBHadron.clone()
matchGenBCHadronB.genJets = cms.InputTag("ak5GenJetsPlusBCHadron", "", "")
matchGenBCHadronC = matchGenCHadron.clone()
matchGenBCHadronC.genJets = cms.InputTag("ak5GenJetsPlusBCHadron", "", "")
GenBCHadronMatchingSequence = cms.Sequence( genParticlesForJetsPlusBCHadron * ak5GenJetsPlusBCHadron * matchGenBCHadronB * matchGenBCHadronC )




# Configuration for matching B-hadrons (jets without neutrinos) ========================================
genParticlesForJetsNoNuPlusBHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusBHadron.injectHadronFlavours = cms.vint32(5)
ak5GenJetsNoNuPlusBHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusBHadron.src = "genParticlesForJetsNoNuPlusBHadron"
matchGenBHadronNoNu = matchGenBHadron.clone()
matchGenBHadronNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusBHadron", "", "")
GenBHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusBHadron * ak5GenJetsNoNuPlusBHadron * matchGenBHadronNoNu )



# Configuration for matching C-hadrons (jets without neutrinos) ========================================
genParticlesForJetsNoNuPlusCHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusCHadron.injectHadronFlavours = cms.vint32(4)
ak5GenJetsNoNuPlusCHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusCHadron.src = "genParticlesForJetsNoNuPlusCHadron"
matchGenCHadronNoNu = matchGenCHadron.clone()
matchGenCHadronNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusCHadron", "", "")
GenCHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusCHadron * ak5GenJetsNoNuPlusCHadron * matchGenCHadronNoNu )



# Configuration for matching B- and C-hadrons (jets without neutrinos) =================================
genParticlesForJetsNoNuPlusBCHadron = genParticlesForJetsNoNuPlusNoHadron.clone()
genParticlesForJetsNoNuPlusBCHadron.injectHadronFlavours = cms.vint32(5, 4)
ak5GenJetsNoNuPlusBCHadron = ak5GenJets.clone()
ak5GenJetsNoNuPlusBCHadron.src = "genParticlesForJetsNoNuPlusBCHadron"
matchGenBCHadronBNoNu = matchGenBHadron.clone()
matchGenBCHadronBNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
matchGenBCHadronCNoNu = matchGenCHadron.clone()
matchGenBCHadronCNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
GenBCHadronMatchingNoNuSequence = cms.Sequence( genParticlesForJetsNoNuPlusBCHadron * ak5GenJetsNoNuPlusBCHadron * matchGenBCHadronBNoNu * matchGenBCHadronCNoNu )


