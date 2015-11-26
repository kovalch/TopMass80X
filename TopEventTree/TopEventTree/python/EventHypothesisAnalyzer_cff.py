import FWCore.ParameterSet.Config as cms

#
# make simple analysis plots for a comparison
# between a simple algorithmic, a gen match and
# an MVA discriminator based event hypothesis
#

# initialize analyzers
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import *
analyzeGenMatch      = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypGenMatch:Key")
analyzeMVADisc       = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypMVADisc:Key")
analyzeHitFit        = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key")
analyzeKinFit        = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypKinFit:Key")

# define sequence
analyzeHypotheses = cms.Sequence(#analyzeGenMatch
                                 #*analyzeMVADisc
                                 analyzeHitFit
                                 )
