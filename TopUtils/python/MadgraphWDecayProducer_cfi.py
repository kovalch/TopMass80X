import FWCore.ParameterSet.Config as cms

MadgraphWDecayProducer = cms.EDProducer('MadgraphWDecayProducer',
    ## Input particle collection of type edm::View<reco::GenParticle>
    src = cms.InputTag("genParticles"),
)
