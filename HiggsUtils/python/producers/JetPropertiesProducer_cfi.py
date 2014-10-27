import FWCore.ParameterSet.Config as cms

jetProperties = cms.EDProducer('JetPropertiesProducer',
	# input jet collection
    src = cms.InputTag("selectedPatJets"),
    svComputer = cms.InputTag("combinedSecondaryVertex"),
)
