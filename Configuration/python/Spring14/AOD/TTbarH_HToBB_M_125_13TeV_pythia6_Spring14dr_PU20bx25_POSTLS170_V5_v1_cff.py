import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/06BD29E8-29D0-E311-85CF-1CC1DE040FA0.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/06E4F119-1DD0-E311-9396-00266CFFC89C.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/0C09E740-17D0-E311-80C6-00266CFFBED8.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/0EF642B6-1AD0-E311-99E4-AC162DA87230.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/1CAB7E58-0BD0-E311-B688-00266CFFBC3C.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/24978213-21D0-E311-84C2-00266CFFCC50.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/30E46872-19D0-E311-AB02-0017A4771010.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/3A972BE1-16D0-E311-A4B9-1CC1DE05B0C8.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/4C8AB75B-15D0-E311-A034-1CC1DE041FD8.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/5438127B-FECF-E311-8A55-D8D385FF0B6A.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/6863A650-32D0-E311-BDFD-0017A4770C3C.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/9ADE8CDE-1FD0-E311-9295-1CC1DE05B0C8.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/A8495261-40D0-E311-9CCB-00266CFFC89C.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/AEDEB939-10D0-E311-BD3D-00266CFFB7D0.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/BE9F481F-25D0-E311-AF29-1CC1DE047F98.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/C21B334F-49D0-E311-9D18-1CC1DE051060.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/CEC98FD1-2DD0-E311-824F-00266CFFC43C.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E8A9688C-18D0-E311-8B39-1CC1DE0437C8.root',
       '/store/mc/Spring14dr/TTbarH_HToBB_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/FCD799F0-4AD0-E311-B2E5-1CC1DE051080.root' ] );


secFiles.extend( [
               ] )

