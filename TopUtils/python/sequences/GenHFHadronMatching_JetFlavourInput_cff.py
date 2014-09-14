import FWCore.ParameterSet.Config as cms


from TopAnalysis.TopUtils.GenHFHadronMatcher_JetFlavourInput_cfi import matchGenHFHadron



# Configuration for matching B-hadrons ================================================================
matchGenBHadron = matchGenHFHadron.clone()
matchGenBHadron.flavour = 5
genBHadronMatchingSequence = cms.Sequence( matchGenBHadron )


# Configuration for matching C-hadrons =================================================================
matchGenCHadron = matchGenHFHadron.clone()
matchGenCHadron.flavour = 4
genCHadronMatchingSequence = cms.Sequence( matchGenCHadron )


# Configuration for matching B- and C-hadrons ==========================================================
matchGenBCHadronB = matchGenBHadron.clone()
matchGenBCHadronC = matchGenCHadron.clone()
genBCHadronMatchingSequence = cms.Sequence( matchGenBCHadronB * matchGenBCHadronC )


# Configuration for matching B-hadrons (jets without neutrinos) ========================================
matchGenBHadronNoNu = matchGenBHadron.clone()
genBHadronMatchingNoNuSequence = cms.Sequence( matchGenBHadronNoNu )


# Configuration for matching C-hadrons (jets without neutrinos) ========================================
matchGenCHadronNoNu = matchGenCHadron.clone()
genCHadronMatchingNoNuSequence = cms.Sequence( matchGenCHadronNoNu )


# Configuration for matching B- and C-hadrons (jets without neutrinos) =================================
matchGenBCHadronBNoNu = matchGenBHadron.clone()
matchGenBCHadronBNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
matchGenBCHadronCNoNu = matchGenCHadron.clone()
matchGenBCHadronCNoNu.genJets = cms.InputTag("ak5GenJetsNoNuPlusBCHadron", "", "")
genBCHadronMatchingNoNuSequence = cms.Sequence( matchGenBCHadronBNoNu * matchGenBCHadronCNoNu )


