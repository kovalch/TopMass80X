import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/08A15648-C110-E411-9CF3-002590A831CA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/0E3D08A9-C610-E411-A862-0025B3E0657E.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/0E6E3678-C510-E411-800F-002590A3C970.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/16E46E35-C010-E411-A0B4-002590A37116.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/1872BE84-4511-E411-B4CD-001E67397D73.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/1CE7D33B-AB10-E411-8BB4-002590A887F0.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/26AD8D51-7317-E411-9229-0025B31E3CC0.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/2CC176FF-C810-E411-9639-001E67396FD1.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/3C152B5F-D211-E411-9C68-001E67396A09.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/3E90250D-8B10-E411-9DE4-002590A371D4.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/4866BDEC-8013-E411-BAC3-002590A80D9C.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/683201DA-FC10-E411-85A6-002590A3C95E.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/6E0E1E31-C010-E411-8F34-002590A37114.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/7A6EA953-9C10-E411-B097-002590A887AC.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/8823E174-7610-E411-9F3E-002590200984.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/8CC5B59B-BB10-E411-9F4B-001E67398B29.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/9099E11B-FA11-E411-AB81-002481E154EC.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/90CD88C6-CE11-E411-8D7C-002590A80DEA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/96AF279A-3D11-E411-BAC0-002590200B60.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/9A3C2D7D-DD10-E411-9AB4-001E67398CA0.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/9EEC30DB-B310-E411-8496-0025B3E05C32.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/A0BF13AD-B610-E411-A49D-F04DA23BBCCA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/A44285CD-AC10-E411-8E2F-001E673989DF.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/B0393652-9210-E411-A13B-0025902009B4.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/B25F5B1C-A810-E411-8CBB-0025901248FA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/BEA05D5F-D811-E411-B625-002590A831CA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/C24B41AC-CC10-E411-975E-002590A371CA.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/C83C48A5-CB10-E411-98A6-001E67396FD1.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/CC36FE9D-1812-E411-8C80-002590A80DD4.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/CC634161-F510-E411-AB14-001E6739665D.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/D02CF7A1-9110-E411-BC4D-002481E75CDE.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/D8DC0E04-1911-E411-B2BB-0025B3E0657E.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/D8FC1042-A510-E411-91EB-002590200814.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E0B3A5C2-A310-E411-A90D-002590A50046.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E43DBBA7-B610-E411-B0FE-002590A37118.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E6394EB5-7012-E411-A97C-002590A88830.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E895F407-E410-E411-8164-0025B31E3C58.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/EEF464C5-B312-E411-B56E-002590200B00.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/F46F8360-D010-E411-89CE-001E67397D00.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/F6936443-0E11-E411-9955-002590A3711E.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/F6FC74FB-0119-E411-B06B-001E67398052.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/FC5ADB67-E710-E411-8D07-001E67398958.root',
       '/store/mc/Spring14dr/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/FE4C49A1-9710-E411-8083-002590200ABC.root' ] );


secFiles.extend( [
               ] )

