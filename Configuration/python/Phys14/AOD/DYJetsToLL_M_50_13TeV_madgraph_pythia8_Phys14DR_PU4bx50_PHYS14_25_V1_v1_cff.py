import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0071F84D-856E-E411-88DA-003048F0E51A.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02224EE7-796E-E411-99C9-00266CF25320.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/04B021D0-826E-E411-9397-0025901D4764.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/08AD2E62-8D6E-E411-84E9-0025904B1428.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0EC8E32F-766E-E411-A913-00266CFFA184.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/104511FD-796E-E411-9AE2-003048D4DEA6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/16582FDD-8F6E-E411-98E2-0025901D4D20.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/242246D9-866E-E411-81AF-003048D43724.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/249C6CE1-796E-E411-8D0C-00266CF32930.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/24A57AB4-806E-E411-9448-002590AC4CC2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2A1BE3CA-806E-E411-ACF4-002590AC4C08.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2EFD9447-766E-E411-80F8-00266CFFA2EC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/340060E6-806E-E411-B9FB-0025904B1428.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/38C98ECC-806E-E411-9B39-003048F0E822.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3C01FDD4-826E-E411-A7C2-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/428E1329-766E-E411-AFA1-00266CFF9FFC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/44E05D49-766E-E411-90F0-00266CFFA184.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A4E1E6A-856E-E411-9741-0025901D4AF4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5A9560C1-856E-E411-B7DF-0025901D4AF4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/602F536D-836E-E411-8A58-002481E0EA06.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/607A6EB4-806E-E411-8C1F-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/608DCBCB-806E-E411-B0E2-0025901D4C92.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/66D7A03D-856E-E411-8636-0025907DC9C6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6A319F4C-886E-E411-9D5D-0025901D484C.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6A6C9DB6-806E-E411-9E44-002590AC4FEC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6AA464BC-806E-E411-8B20-00266CF3338C.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6E2FE7B4-806E-E411-8CB6-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/72933E06-896E-E411-B825-0025901D4764.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7498F028-766E-E411-837D-003048F0E55A.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/76EFF56E-816E-E411-B045-003048F0E186.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7A72B932-A06E-E411-847B-003048F0EBBE.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7A7CCB45-766E-E411-A7B5-00266CFFA704.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7C2E28A7-816E-E411-957A-002590DB9252.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/804DB921-766E-E411-98C4-00266CF422D8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/82EE62C7-796E-E411-A0AA-008CFA104E64.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8680EAC5-806E-E411-98D0-00266CFFA5C4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8852BAA8-816E-E411-88E5-0025907DCA72.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8AFE17DF-8A6E-E411-85C7-0025901D4B20.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8EEA8C5E-886E-E411-B320-002481E15176.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/90072AF7-796E-E411-BF32-002590AC4CA6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A2873F0D-716E-E411-8582-00215AD4D6AE.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AA4018ED-796E-E411-AC6C-00266CF2AACC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AC0F022C-766E-E411-9EAD-0025907FD430.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AE538E59-856E-E411-A3C9-002590AC4FEC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AE6D6A6B-6E6E-E411-B883-002590DBDFE0.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B82155DD-806E-E411-B0ED-002590AC4C26.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BA1093BD-806E-E411-B428-002481E945C4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BA7C6D01-926E-E411-A474-00266CF330D8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BEBAD5DF-866E-E411-84B1-003048F0E194.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C04742B8-806E-E411-A3CA-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C06CEAC8-706E-E411-B715-0025907DBA06.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C6ADDA55-B66E-E411-B92E-002590AC501A.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CA9467C3-806E-E411-BCA1-00266CF33054.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CC8AC198-8B6E-E411-B654-003048F02CB6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D8119980-B66E-E411-AAC2-0025901D4B02.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E09182EB-896E-E411-845D-002590DB91D2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E208E131-766E-E411-A0AE-003048F0E3B0.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EAFBE7B6-806E-E411-B2EB-002590DB927A.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F41467C6-796E-E411-8042-00266CFFA418.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F680E946-866E-E411-86F4-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FE0428C5-886E-E411-8E81-002481E0D69C.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FE9ACF2B-766E-E411-AEC1-002590AC4CA6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/0473A12E-736E-E411-B0B4-002590AC4BF6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/04D19299-B86E-E411-AB5F-003048D436D2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/04D54243-8D6E-E411-9FB7-0025904B12E2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/0ABDA828-7C6E-E411-A199-003048D2BB22.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/10564669-816E-E411-8D5B-002590AC4C48.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/12B7EB4A-976E-E411-9FAB-00215AEDFC8E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/12BE2FBA-736E-E411-BF94-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/162A980B-7C6E-E411-9E3C-002590AC5062.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/189A673E-856E-E411-9F2A-002590AC4CC4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/1AD6A62A-766E-E411-8A17-00266CFFA5E0.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/2448CB98-856E-E411-9A55-003048D46122.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/24595206-8A6E-E411-B5AF-0025907DC9D6.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/2AE1AF8D-716E-E411-AD71-0025907FD2DC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/2EBFC25B-816E-E411-902C-003048F0E3AE.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/2ECE054B-816E-E411-93AA-002590DB9296.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/32BF1682-8A6E-E411-86FF-002590AC4C74.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/34AEF709-906E-E411-A27B-002590AC539E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/3800795F-816E-E411-8183-003048D437BA.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/38C1F72F-856E-E411-8439-002590D94F8E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/48006512-876E-E411-86DF-003048F0E3B4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/4A0568F3-966E-E411-9358-003048F0E780.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/540FDC0E-856E-E411-9D01-002481E100F8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/54DEE875-6F6E-E411-9AD3-0025904B12E2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/5A9793CB-786E-E411-A103-002590494C98.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/5CC431EA-786E-E411-AA57-003048F0E5A2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/5E8D9369-816E-E411-8E15-002590AC4CB8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/64143C5A-746E-E411-A627-002481E7628E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/72DABCCE-786E-E411-B92D-0025904B5FB8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/76F9E9B6-B86E-E411-8B1E-002590AC4C7C.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/8856C66B-6D6E-E411-9380-0025907DC9AC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/8A018B59-816E-E411-9F6E-002590AC4C24.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/904C228B-716E-E411-B34C-002590DB918A.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/9C5494F1-846E-E411-A4ED-002590AC5082.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/A60748FE-706E-E411-A218-0025907FD242.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/ACEE2D6D-896E-E411-983B-002590AC5082.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/AE4BDAA3-886E-E411-A5A6-002590D94F8E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B01C954D-746E-E411-9941-0025901D40CA.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B2AEF2A8-866E-E411-A658-0025904B0FE4.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B48C0EDC-846E-E411-A0F2-002590DB9182.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B4EA767A-756E-E411-AD8E-00266CFFA2B8.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B819C2FD-876E-E411-ABD4-002590AC4C74.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B8944DEB-846E-E411-9040-00266CFFA7C0.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/BA511F45-816E-E411-A14E-002590AC5482.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/BA6ED8E8-846E-E411-95AD-003048D437BA.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/BAFA5D33-6F6E-E411-83B6-0025907FD442.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/C0A229E5-846E-E411-8622-0025904B1026.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/C6BDE903-856E-E411-B43D-003048F0EBBE.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/DEB29803-886E-E411-8FD6-002590AC4C8C.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/E40CD2E5-846E-E411-B35F-0025901D4C92.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/E8F8D8A0-6E6E-E411-B96C-0025907FD3CC.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/ECA8A725-856E-E411-9278-D8D385FF6C5E.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/F0D0770F-7C6E-E411-9BC4-00266CFFA780.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/F4D218E7-846E-E411-A1C5-003048D47912.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/F6673B3B-866E-E411-B049-002590DB91D2.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/F67A9437-8C6E-E411-8111-002590AC4D32.root',
       
'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/FAAF160F-7C6E-E411-8450-00266CF3322C.root' 
] );


secFiles.extend( [
               ] )

