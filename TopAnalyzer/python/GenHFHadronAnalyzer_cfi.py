import FWCore.ParameterSet.Config as cms

analyzeGenHFHadron = cms.EDAnalyzer("GenHFHadronAnalyzer",
    genJets = cms.InputTag('ak5GenJets','','SIM'),   
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(False),
    noBBbarResonances = cms.bool(True),
    doValidationPlotsForImprovedHadronMatching = cms.bool(True)
)



