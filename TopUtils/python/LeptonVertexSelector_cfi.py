import FWCore.ParameterSet.Config as cms

leptonVertexSelector = cms.EDProducer('LeptonVertexSelector',
    
    # Input sources for electrons, muons and vertices
    electrons = cms.InputTag("selectedPatElectrons"),
    muons = cms.InputTag("selectedPatMuons"),
    vertices = cms.InputTag('goodOfflinePrimaryVertices'),
    
    # Selection cuts on dxy and dz, with respect to the leading primary vertex
    # If value is negative, no selection is applied for this
    electronDxyMax = cms.double(-999.),
    electronDzMax = cms.double(-999.),
    muonDxyMax = cms.double(-999.),
    muonDzMax = cms.double(-999.),
)
