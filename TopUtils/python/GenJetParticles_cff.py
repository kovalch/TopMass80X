import FWCore.ParameterSet.Config as cms
from TopAnalysis.TopUtils.GenJetParticles_cfi import *


genParticlesForJetsNoNuPlusNoHadron = genParticlesForJetsPlusNoHadron.clone()
genParticlesForJetsNoNuPlusNoHadron.ignoreParticleIDs += cms.vuint32( 12,14,16)