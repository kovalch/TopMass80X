import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/0E5DBBD5-606C-E411-8E1F-002590A37116.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/16816544-586C-E411-BF87-002590A80DDA.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/1A67568C-586C-E411-8C45-002590A3A3D2.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/2C96D02B-586C-E411-B52A-001E6739811A.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/304CEB95-6E6C-E411-9A16-002481E150C2.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/32B6450F-616C-E411-81F2-002590A887EE.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/42A6D4BA-586C-E411-A52A-002590A36FB2.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/5673433E-586C-E411-B12A-002590A8881E.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/582363B1-856D-E411-92CE-002590200874.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/5CE6F51E-5B6C-E411-AC4E-002590A80E1E.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/665EC632-596C-E411-9AF7-002590A80D9C.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/7A2A56F7-606C-E411-92FF-002590A37122.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/84CAD93A-596C-E411-A129-002590A8881E.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/8848D8CC-6E6C-E411-90D0-002590A371AA.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/928359D1-586C-E411-9031-002590A83354.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/9A03D73F-586C-E411-9A5C-0025902009E4.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/ACB4D3A7-586C-E411-BA72-002590A36FB2.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/BC6AB861-766C-E411-97E3-002590A88806.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/BEDF7E17-5B6C-E411-AAA9-002590A8881E.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/C4986348-5E6C-E411-933D-002590A3C954.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/CE079B24-596C-E411-9FAA-001E6739811A.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/D0E7AB34-606C-E411-9B9A-002590A80DFA.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/E206C4E3-5D6C-E411-A02A-001E67398CE1.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/00000/E8A1F176-606C-E411-8C30-001E6739670C.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/08C9EF94-846C-E411-86A1-002590A3C984.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/0A7DC29C-856C-E411-9638-0025B3E06578.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/105D8FF8-486D-E411-AD7E-002590A8882C.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/20BB74B9-856C-E411-8F17-002590200B38.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/2436EA45-6E6C-E411-8B18-001E673982E6.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/301E5A01-496D-E411-AC0E-002590A3C96C.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/3CCE3237-756C-E411-92FB-002590A831CC.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/40B9A1B0-846C-E411-BB4A-0025B3E063BA.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/4EF16431-6E6C-E411-BFEB-002590200894.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/6862209A-6E6C-E411-9B28-001E67397233.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/709741E7-6D6C-E411-9738-002590200B68.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/7AF771AD-856C-E411-B075-002590A831B4.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/7C1566D2-836C-E411-9D3D-001E67398BE7.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/84DAAF67-6E6C-E411-B82B-002590A4FFE8.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/8C64F68A-6D6C-E411-A152-002590A3C978.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/90355D78-6C6C-E411-A515-002590A4FFE8.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/9CF1D9D7-6D6C-E411-8875-D8D385FF7678.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/B8A43275-6C6C-E411-8B28-D8D385FF7678.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/C218B561-756C-E411-A144-002590A3C978.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/FC71A40F-6D6C-E411-AE71-002590A4FFE8.root',
       '/store/mc/Phys14DR/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v1/10000/FCC6FF4A-6D6C-E411-8A8C-001E673982E6.root' ] );


secFiles.extend( [
               ] )

