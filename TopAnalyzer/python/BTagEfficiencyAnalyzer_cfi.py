import FWCore.ParameterSet.Config as cms

analyzeBTagEfficiency = cms.EDAnalyzer("BTagEfficiencyAnalyzer",
  ## jet collection (after jet selection, before b-tagging)
    jets  = cms.InputTag("tightLeadingPFJets"),
  ## vertex collection                                     
    vertices = cms.InputTag("goodOfflinePrimaryVertices"),
  ## name of b tag algorithm in pat::Jet
    bTagAlgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
  ## b tag discriminator cut value (here: medium working point)
    bTagDiscrCut = cms.double(0.8),
  ## set to 1 if terminal text output is desired
    #verbose  = cms.int32(  0)
  ## binning for b and c jets:
  ## bin edges for variable bins: vector of low-edges for each bin (including overflow), i.e. nBins+1 entries
    binsPtB     = cms.vdouble(  0., 10., 20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 670., 1000., 2000.), ## binning of BTV
    binsEtaB    = cms.vdouble(  0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3. ),
  ## binning for light jets
  ## bin edges for variable bins: vector of low-edges for each bin (including overflow), i.e. nBins+1 entries
    binsPtL     = cms.vdouble(  0., 10., 20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 670., 1000., 2000.),
    binsEtaL    = cms.vdouble(  0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3. )
)