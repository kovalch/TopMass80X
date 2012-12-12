import FWCore.ParameterSet.Config as cms

#
# make simple analysis plots for a comparison
# between a simple algorithmic, a gen match and
# an MVA discriminator based event hypothesis
#

# initialize tree creator
createEventTree = cms.EDAnalyzer("EventTreeCreator")

# initialize analyzers
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import *
analyzeGenMatch      = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypGenMatch:Key")
analyzeMVADisc       = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypMVADisc:Key")
analyzeHitFit        = analyzeHypothesis.clone(hypoClassKey = "ttSemiLepHypHitFit:Key")

# define sequence
analyzeHypotheses = cms.Sequence(createEventTree
                                 #*analyzeGenMatch
                                 #*analyzeMVADisc
                                 *analyzeHitFit
                                 )
