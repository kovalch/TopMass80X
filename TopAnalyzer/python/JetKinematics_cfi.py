import FWCore.ParameterSet.Config as cms

analyzeJetKinematics = cms.EDAnalyzer("JetKinematicsAnalyzer",
    ## input collection                             
    src = cms.InputTag("selectedLayer1Jets"),
    ## event weight
    weight = cms.InputTag("eventWeight"),
    ## use the weight or not                             
    useWeight = cms.bool(True),
    ## special parameters for muon quality analysis
    analyze   = cms.PSet(
      ## fill histograms for 1.,2.,3.,... leading
      ## jet? -1 corresponds to 'all'
      index = cms.int32(-1),
      ## correct jets up to the correction level given
      ## at SWGuidePATDataFormats#pat_JetCorrFactors;
      ## the expected form is of type 'step' or 'step:
      ## flavor'
      correctionLevel = cms.string("abs")
    )                                       
)



