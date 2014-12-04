import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/067E448A-4270-E411-A1E1-002618943959.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/06FD7192-4270-E411-ACDA-0025905A6070.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0E322D90-4270-E411-9EF3-003048FFCB84.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/10762291-4270-E411-AB3E-0025905B85A2.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1E943390-4270-E411-90E6-0025905B860E.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/46744990-4270-E411-AACC-0025905A6090.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4AFE438A-4270-E411-A98F-003048FFD754.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/521B8290-4270-E411-8C82-0025905A6076.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/82761F90-4270-E411-B910-0025905A60B6.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9EF8CF8C-4270-E411-8A1C-0025905A60F4.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A6EB748D-4270-E411-8078-0025905A48FC.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AE8F6A8C-4270-E411-A041-002618943913.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C2DB8F8E-4270-E411-94FF-0025905A6082.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CC902A8A-4270-E411-8AA1-0025905A6056.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CE547F8D-4270-E411-9D33-0025905A48D6.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CEA1718A-4270-E411-B68A-0025905A613C.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CEB10D8E-4270-E411-9EC2-0025905AA9CC.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CED7FE8F-4270-E411-B2FC-003048FFD75C.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F65C908D-4270-E411-9499-0025905B8572.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FEBE838F-4270-E411-A003-0025905938D4.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/10000/32885C77-1B70-E411-9711-0025905B85AA.root',
       '/store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/10000/DEEEA278-1B70-E411-8F41-00261894389C.root' ] );


secFiles.extend( [
               ] )

