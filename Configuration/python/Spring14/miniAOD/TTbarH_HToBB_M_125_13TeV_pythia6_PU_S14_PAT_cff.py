import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       'file:/nfs/dust/cms/group/topcmsdesy/MiniAOD/Spring14/DataSets/TTbarH_HToBB_M-125_13TeV_pythia6_PU_S14_PAT.root' ] );


secFiles.extend( [
               ] )

