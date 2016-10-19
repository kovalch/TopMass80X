import FWCore.ParameterSet.Config as cms

####################################################################################################################################################
# Dataset: /TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/MINIAODSIM      #                                         
# Events : 182123200                                                                                                                                  #
# xSec   : 831.76pb (NNLO)                                                                                                                         #
# eff    : 1.0                                                                                                                                     #
####################################################################################################################################################

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring( *(
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/0030C76E-1F29-E611-BF34-0CC47A4D763C.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/3CBFA91F-C128-E611-A4D9-0025905A607A.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/3CD094E1-DA28-E611-9785-0022195000BF.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/3E018595-F928-E611-B54C-002481CFE834.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/3E842C17-C628-E611-9D17-008CFA1979B0.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/7C354BBE-FC28-E611-A631-6C3BE5B533A8.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/7CC5950A-F928-E611-B609-ECB1D79E5C40.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/7CDA155A-0029-E611-AF4B-0025905C42A8.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/BC64B06F-C228-E611-9EFD-0025905B85B2.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/BC9BA042-BE28-E611-A745-A0000420FE80.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/FCB5C8C7-E828-E611-B2BA-002481DE4BF0.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00001/FE0C4673-DF28-E611-B883-848F69FD46BB.root'
    ) )
)
