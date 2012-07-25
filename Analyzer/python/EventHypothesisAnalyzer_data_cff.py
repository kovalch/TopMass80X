import FWCore.ParameterSet.Config as cms

#
# make simple analysis plots for a comparison
# between a simple algorithmic, a gen match and
# an MVA discriminator based event hypothesis
#

# initialize analyzers
from TopMass.Analyzer.EventHypothesisAnalyzer_cfi import *
analyzeMVADisc       = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypMVADisc:Key", mets = "patMETsPF")
analyzeHitFit        = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key", mets = "patMETsPF")

# define sequence
analyzeHypotheses = cms.Sequence(analyzeMVADisc *
                                 analyzeHitFit)
