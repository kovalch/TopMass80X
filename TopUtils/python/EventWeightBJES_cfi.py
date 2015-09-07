import FWCore.ParameterSet.Config as cms

EventWeightBJES = cms.EDProducer("EventWeightBJES",
  genJets               = cms.InputTag("ak4GenJets"),
  genParticles          = cms.InputTag("genParticles"),
  nuDecayFractionSource = cms.double(0.25),
    ## Pythia6 is 0.25
  nuDecayFractionTarget = cms.double(0.25),
    ## use up:0.262, down:0.238, based on PDG
  fragSourceFile        = cms.FileInPath("TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star.root"),
  fragTargetFile        = cms.FileInPath("TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star.root")
    ## for variation: TopAnalysis/TopUtils/data/MC_BJES_TuneZ2star_rbLEP.root
)
