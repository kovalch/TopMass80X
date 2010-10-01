import FWCore.ParameterSet.Config as cms

filterMuonVertexDistance = cms.EDFilter("MuonVertexFilter",
    vertices = cms.InputTag('offlinePrimaryVertices'),
    muons    = cms.InputTag('selectedPatMuons'),
    
    ndof   = cms.uint32(4),
    rho    = cms.double(2),
    z      = cms.double(24),
    z_dist = cms.double(1)   
)
