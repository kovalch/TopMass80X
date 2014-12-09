import FWCore.ParameterSet.Config as cms

jetProperties = cms.EDProducer('JetPropertiesProducer',
	# input jet collection
    src = cms.InputTag("selectedPatJets"),
    pfCands = cms.InputTag("packedPFCandidates"),
    svComputer = cms.InputTag("combinedSecondaryVertex"),
)
