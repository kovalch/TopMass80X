import FWCore.ParameterSet.Config as cms

matchGenHFHadronLocal = cms.EDProducer("GenHFHadronMatcherLocal",
    genJets = cms.InputTag('ak5GenJets','','SIM'),   
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(False),
    noBBbarResonances = cms.bool(True),
)



