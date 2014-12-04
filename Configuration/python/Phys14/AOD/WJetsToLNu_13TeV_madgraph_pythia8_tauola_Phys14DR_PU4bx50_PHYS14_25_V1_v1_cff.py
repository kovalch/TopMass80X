import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/00DBC651-E86F-E411-88A4-0026189438FF.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/022DE69B-CB6F-E411-9BCC-00261894385A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02693744-CD6F-E411-A8EB-002590593878.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02DDA1A5-CE6F-E411-91EB-0025905A60B4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/041313AA-CB6F-E411-9BB1-00248C0BE01E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/04B7C83E-D86F-E411-89A0-0026189438C9.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/060286FB-E36F-E411-9417-003048FFD7C2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/061F4BBA-CF6F-E411-ACF3-002618943939.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0657C92D-D36F-E411-80C9-0025905B85A2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/06FD6BB0-CB6F-E411-8F65-0025905A48BC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0808EA1C-CF6F-E411-8BD8-0025905A48D0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/081FEC5F-CD6F-E411-AC2E-0025905A6082.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/084E7EF6-B66F-E411-9648-0025905A6064.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0A50ADB9-CF6F-E411-857F-0025905964BE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0ADACEF9-DE6F-E411-8299-002618943896.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0CBE8F9D-CB6F-E411-A507-0025905A612A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0CD46377-CB6F-E411-B735-0025905B8572.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0E442E28-D46F-E411-9619-0025905AA9CC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0EC6B640-D96F-E411-B131-003048FFD75C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/14504D0B-CF6F-E411-ACE2-002618943951.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/14970FF6-D06F-E411-B569-002590593878.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/160CFBC6-CA6F-E411-A9C0-002618943856.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/18515D52-CD6F-E411-B19A-0025905A48C0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/18DD2A65-B16F-E411-AE88-002618943849.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1A157299-B56F-E411-8ED6-0025905A6122.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1A475751-CD6F-E411-8F8F-0025905A48B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1AC48B7E-3970-E411-B601-0025905A48D6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1AE67B42-CD6F-E411-A9F3-0025905A48D0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1C4287A3-CB6F-E411-A637-003048FFD76E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/1E98B7F8-DE6F-E411-9D58-002590593902.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/22321304-D06F-E411-AABF-0025905A60D6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/22E6A9E8-D26F-E411-8E64-003048FFCB9E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/24AAF8A8-CB6F-E411-9E64-00261894393D.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/24F33EBC-FB6F-E411-9232-00261894397E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/26975C8D-AD6F-E411-94B6-0025905A6132.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/26C6540D-CF6F-E411-8B1F-002618943947.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/280B7665-CD6F-E411-AF2D-003048FFCB6A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/285FBEC6-CF6F-E411-8EDD-0025905A60B4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2AE50C5F-DE6F-E411-A80F-002618943949.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2C561A46-EF6F-E411-B857-0025905B8592.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2C97971B-CF6F-E411-9EA7-0025905A6068.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2E5D9894-B56F-E411-BA25-003048FFD7C2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/2EF95CAF-CB6F-E411-B427-002590593876.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3047A97C-AD6F-E411-BBAA-0025905A60BC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/30685871-DE6F-E411-B193-0025905B8562.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3433896C-DE6F-E411-A828-002590593920.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/380E15AB-CB6F-E411-BED8-002590593878.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/38943E66-B16F-E411-AA26-0025905A6088.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3C073CB5-CF6F-E411-AC02-0025905A609E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3C96D415-CF6F-E411-ABB2-0025905A48FC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3CC12BB1-CF6F-E411-A16F-002590593902.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3E1BD460-CB6F-E411-B569-002590596490.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/3EF5A4C5-CF6F-E411-A445-002590596498.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4032A3A0-CB6F-E411-A9CA-003048FF86CA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/408C2BDB-D16F-E411-BBFD-002618943933.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/40F6099B-B56F-E411-9119-003048FFCB8C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/42A83920-CF6F-E411-91FE-003048FFD75C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/42AFCBF3-D06F-E411-A650-0026189438A0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4479F453-B36F-E411-AA4F-0025905A610C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/44C3A563-CB6F-E411-BAFC-0025905B85D0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/46644FB8-CF6F-E411-BEC4-0025905A48EC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/48F65356-D06F-E411-8F4F-003048FFD7BE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A5556A7-CB6F-E411-AFAB-0026189438D8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A5671FA-E36F-E411-BCED-0025905B8576.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A6D294C-D76F-E411-B873-0025905A608E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A71174D-CB6F-E411-8EF8-003048FFCC1E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A7D9841-D16F-E411-A707-0025905A60D2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4A9159B9-CF6F-E411-8F70-0025905964A6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4AA747FB-E36F-E411-A052-003048FFD770.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4ACFE0C6-CF6F-E411-8216-0025905A60B4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4CB1E0BB-FB6F-E411-89C8-0026189437FD.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4CF9663E-D86F-E411-AD5B-0025905B860E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4E03914B-D16F-E411-A6A5-0025905A60EE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4E6DF084-C96F-E411-8936-0025905A4964.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/4E8BEC0A-CF6F-E411-98AE-002618943970.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/506121E3-D06F-E411-81FD-002618943876.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/50EB29A1-CB6F-E411-968A-0025905B855C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/521EF5CE-CF6F-E411-8886-0025905A606A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/52539D15-CF6F-E411-89A8-0025905A60EE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/528471DB-EA6F-E411-AD6F-0025905A611C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5410027D-3970-E411-9E53-0026189438E7.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5446F8CE-CF6F-E411-9ADC-0025905A606A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5456CB83-C96F-E411-B958-0025905AA9CC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/56146801-D16F-E411-8DDC-0025905A60CE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5661C996-B56F-E411-9C82-0025905A610A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/583A9C44-CD6F-E411-8DEA-00259059642E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/58D8FBFD-DE6F-E411-9FBA-003048FFD796.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/58E280BF-CF6F-E411-8BFA-0025905B85D8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/58E7C4E4-E36F-E411-815F-002618943918.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5A2301C5-CF6F-E411-B14B-0025905A6134.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5A305C9F-B56F-E411-A04D-003048FFD7BE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5A7CDB4D-CD6F-E411-BB7B-003048FFD732.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5C441C79-CE6F-E411-B94D-003048FFD730.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5C7F28C2-CA6F-E411-98B7-002618943867.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/5CC718E4-CA6F-E411-B01E-0025905964C4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/60C1CFBB-CF6F-E411-936D-0025905A48D6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/60EBFD40-AC6F-E411-87DA-0026189438DE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6288E267-DE6F-E411-88A5-0026189438DF.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/64E540D3-B56F-E411-9DF5-0025905A612A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6815307E-C96F-E411-853F-003048FFD754.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/682B4B2E-B66F-E411-A0B3-0025905964BA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/68766047-D86F-E411-8437-0025905B8596.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/689FC0C1-CF6F-E411-9E35-0025905B8576.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6A3CB659-CB6F-E411-B5EA-00261894386B.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6AB2EE79-DE6F-E411-BE9C-0025905A608E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6ABFAC14-D86F-E411-9238-0025905B857C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6ACEF3AC-CB6F-E411-B087-0025905A6134.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6C29B370-CD6F-E411-97DC-0025905B8606.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6CF79551-D76F-E411-9350-003048FFCB84.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/6EEF753E-D96F-E411-A906-0025905B858E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7020E671-CB6F-E411-BFFE-0025905A6118.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/703787D3-CA6F-E411-B790-0025905A60D0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/707AF0CC-CA6F-E411-983E-002590596468.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/708EB75A-CB6F-E411-A2F1-0025905A60CA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/70F4CCA7-CB6F-E411-9DDF-002618943866.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/727AE567-DE6F-E411-A1EF-0025905B860C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/746AC140-D36F-E411-BDA5-0025905A6084.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7656AF2C-CA6F-E411-B0A6-003048FFCB74.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/766CFC69-B16F-E411-A823-0025905A608E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/76B0F596-B56F-E411-BA1C-0025905A609E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/788C195C-CB6F-E411-8F17-0025905A607A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/78E098BB-CF6F-E411-8570-0025905A60B6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/78E432B9-CF6F-E411-A746-002618943985.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7CBDE25C-E56F-E411-96DE-002618943915.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7EE81620-B76F-E411-ABD8-0025905B861C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/7EF2FDE5-D26F-E411-9B29-0025905B861C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/800286FD-D06F-E411-BCF0-0025905A60CE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8206449E-CB6F-E411-B3BB-002618943852.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/82858CD1-AD6F-E411-9146-00261894383B.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/82B3C0DB-D26F-E411-A128-0025905A6122.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/82CD376E-CD6F-E411-A6A4-0025905B860E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/845107C0-CB6F-E411-8251-0025905A60DA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/84D03BF9-D46F-E411-BB07-0025905B855E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/84FDACFA-D06F-E411-8826-0025905A60B8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8678AD3B-D16F-E411-84EE-003048FFCB9E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/867F1FA5-CB6F-E411-B21D-0025905964C0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/86F940F6-DE6F-E411-BB21-0025905A612A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8835DFAF-CB6F-E411-A496-0025905AA9CC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/885F11DC-D46F-E411-B5A7-0025905A48E4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8A0B782C-CA6F-E411-8AD0-0025905A60D2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8A4BF652-B36F-E411-AF9C-0025905938A8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8C1A9C34-F16F-E411-93EA-0025905A60D2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8C579E9F-CE6F-E411-A5B3-0025905A60EE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8C887BFD-D06F-E411-BDF5-002618943856.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8E087BB9-D46F-E411-9434-0025905B8598.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/8EE63240-CC6F-E411-831F-0025905A6080.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/900F7362-CB6F-E411-BC6D-0025905A48F0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/90977638-D86F-E411-A6FA-002618FDA210.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/90C49EC3-CF6F-E411-9726-0025905B858A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/90DFB3F7-DE6F-E411-9334-003048FFCBA4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9228E255-CB6F-E411-BE38-003048FFCB96.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/92DC9E6D-DE6F-E411-A864-0025905B8592.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/94141F96-C96F-E411-9174-0025905A610A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/94EE0723-CA6F-E411-9458-002618943867.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/96C3CEC9-CB6F-E411-8906-0025905A6084.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/96C6BB5D-DE6F-E411-A9BA-00261894394D.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9869C009-CF6F-E411-98B4-002618FDA237.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/98ED8631-E46F-E411-96BC-002618FDA250.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9A6B7B7D-CE6F-E411-AB8E-0025905A612C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9A853BCF-CA6F-E411-8451-0025905A60BE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9AD45F11-CF6F-E411-B113-0025905A6082.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9CA76CCA-CA6F-E411-B2A9-0025905A6060.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9CD581D5-D16F-E411-A55D-0026189438E8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9E02BAF7-DE6F-E411-9477-002590593902.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9EE34C43-CA6F-E411-A4AB-0025905A48B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/9EE8F86E-DE6F-E411-9E8F-0025905938AA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A033EAA7-CB6F-E411-A446-002618943849.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A2E8B77A-CE6F-E411-9360-0025905A608A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A44F4B97-B56F-E411-91EA-003048FF9AA6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A64C5578-CB6F-E411-94B9-0025905A60E0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A6B55D57-CB6F-E411-A81A-003048FFD720.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A8A33FDE-EA6F-E411-9BC4-0025905A60B6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/A8C2ECAE-CB6F-E411-A478-003048FFD79C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AA5F8D5C-DE6F-E411-8DEB-00261894398A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AC087D21-CA6F-E411-B060-00261894389D.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AE1F4D7A-C96F-E411-AF2C-0025905A60FE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AE6F0A9D-CB6F-E411-99B3-002618943980.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/AEA57953-B16F-E411-B8A6-002618943970.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B080FA3E-D16F-E411-BD77-00259059642E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B0835578-CB6F-E411-BAD7-0025905A60E0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B083CF82-C96F-E411-9D00-0025905A6110.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B2DC7409-CF6F-E411-B9BA-002618FDA237.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B4854218-D46F-E411-9007-0025905B8576.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B67D6C95-B56F-E411-84B7-0025905A607E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/B6DD5D85-CE6F-E411-9B5E-0025905938A8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BE68D5BF-CF6F-E411-B0BB-0025905964B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BE6C174A-CD6F-E411-AD71-0025905A48BC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BE9BB0E2-D06F-E411-94DF-002618943915.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/BEA6ECD8-CA6F-E411-9546-0025905A48D0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C20605B1-CB6F-E411-995F-0025905A612C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C445FCF1-CF6F-E411-9095-0025905A60E0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C4C91277-C96F-E411-A5E0-002618943845.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C4F9338A-CE6F-E411-9E36-0025905A60B6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C842220D-CF6F-E411-B22A-0025905A60FE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/C8F0153A-D86F-E411-9475-0025905B85B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CA313BDC-D26F-E411-8230-00259059649C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CACEBFD0-AD6F-E411-AC18-002618943919.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CAD6E3D4-CD6F-E411-A283-0025905A48E4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CAF186FB-E36F-E411-9AF3-003048FFCC2C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CC1DF307-CF6F-E411-A2E2-002618943900.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/CEAB013C-AC6F-E411-A74E-00261894393B.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D081BB14-E76F-E411-8818-003048FFD728.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D20252A8-CB6F-E411-9C9B-00261894397F.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D260562E-CF6F-E411-94C6-0025905A607E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D299C289-C96F-E411-AFC7-0025905B859E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D2DB3529-CA6F-E411-8306-0025905964CC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D66B69F2-B66F-E411-8B5B-0025905B8596.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/D6DFDEB3-CF6F-E411-A501-0026189437FE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DA59FA5B-DE6F-E411-B25C-002618943975.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DA9CBBF7-DE6F-E411-87A5-0025905A60B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DE1E50C8-CB6F-E411-9C1E-0025905A60E0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DE3F7227-CA6F-E411-8B36-00259059642E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DEBC2827-DD6F-E411-9A8D-0025905A6076.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/DEDFE4E2-E26F-E411-AA90-0025905A6104.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E0D39121-CF6F-E411-9EA0-0025905A60B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E2EF586E-B16F-E411-9FD9-0025905A48D8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E4956795-B56F-E411-A070-0025905938AA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E62EB642-D86F-E411-A916-00261894388A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E85A3F9D-CB6F-E411-9F60-002618943954.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E863CD0E-D26F-E411-A3E9-0025905A6136.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/E8A69459-CD6F-E411-83C8-0025905A6138.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EA354543-D96F-E411-B9BB-0025905A48BA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EA401286-D26F-E411-9190-0025905B855C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EA6537C1-CF6F-E411-8B21-00259059649C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EAC6EF68-CD6F-E411-AB9A-0025905B858A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EC4F6F16-D46F-E411-8669-002618943975.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/ECC13A6A-DE6F-E411-8FB3-0025905A48BA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/ECD6B647-B66F-E411-94AF-0025905A48C0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EE41F6B2-CB6F-E411-9F47-0025905B860E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EEA6003C-D86F-E411-8B5D-0025905A60BE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/EECB34AE-CB6F-E411-ABDF-0025905938D4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F029EA0D-CF6F-E411-A850-0025905A60C6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F49A723E-D86F-E411-8705-0025905B855C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F4D635A1-CB6F-E411-AA6D-002618943956.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F6276C5E-CB6F-E411-918B-0025905B85AA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F6AC30C0-CF6F-E411-9F86-002618943809.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F874573C-D16F-E411-9698-0025905B8590.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/F8943340-CC6F-E411-A9B7-0025905A6080.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FA2AEF60-CB6F-E411-85EA-0025905A6134.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FCC1F616-CF6F-E411-B7B2-003048FFCC18.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FE08486D-DE6F-E411-B647-003048FFD71E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FE1E8CB3-CF6F-E411-A5E8-00259059649C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FEE87C0D-CF6F-E411-8C09-0025905A60C6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/00000/FEE96FA9-CB6F-E411-B323-0025905A611C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/00AE2020-E36F-E411-9CB2-0025905B85E8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/00E19B31-D36F-E411-9BF7-00261894390E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/021A04B9-DF6F-E411-9A86-002618943924.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/041FCFCB-D76F-E411-AD86-002618943979.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/088FF2D9-F06F-E411-BE3B-003048FFD744.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/1E67F06D-E96F-E411-A9FB-0025905A60F8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/2847BA10-D96F-E411-80A7-0026189438AC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/480ACB51-DB6F-E411-8D60-002618943866.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/52651526-2570-E411-8608-003048FFCBA4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/586D4326-2570-E411-9309-003048FF9AA6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/608F920B-E56F-E411-A43E-002618FDA204.root' 
] );
readFiles.extend( [

       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/6873E4BF-E86F-E411-ADE3-0026189438D7.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/68F52D72-DE6F-E411-8D49-0025905A60E4.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/6A6E5914-D26F-E411-81F1-0025905A60B6.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/824701D4-D36F-E411-B1E6-002618943836.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/8267053E-E16F-E411-8AE0-0025905A48B2.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/8E5922C8-E66F-E411-AE26-00259059649C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/9ACFA5DB-E16F-E411-B158-0025905A60AA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/9C83F0A3-F56F-E411-9967-0025905938A8.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/9EC49FD9-E16F-E411-912A-0025905A613C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/A855E72C-D46F-E411-A11C-003048FFD76E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/AE408FE5-DC6F-E411-B46D-0025905A6122.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/AEFABE1C-DF6F-E411-A647-0025905B861C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B2F147A5-DC6F-E411-813F-0025905A60DE.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B4A96F6A-D66F-E411-B01D-0025905A60AA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/B69D59DD-DD6F-E411-BA82-00261894383B.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/C4766040-D36F-E411-A0F0-0025905A606A.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/C4855DAE-E56F-E411-971F-0025905A607E.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/C8F33C53-DB6F-E411-8139-002618943809.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/CC432113-D96F-E411-8040-002590596490.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/CEBE98E2-DD6F-E411-A7E2-0025905A48BA.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/D0042E31-EE6F-E411-9945-0025905A608C.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/D4F1CB20-E36F-E411-9BBC-0025905B8610.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/D60D8BC4-DF6F-E411-9982-0025905A60B0.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/D8FEAEA2-E06F-E411-9C08-002590593902.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/DC3B91ED-F46F-E411-A460-0025905A6066.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/EE253074-DD6F-E411-9208-0026189438AD.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/F8E9D427-EE6F-E411-9BBD-002354EF3BDC.root',
       
'/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU4bx50_PHYS14_25_V1-v1/10000/FC35C579-DE6F-E411-9A9D-002618943947.root' 
] );


secFiles.extend( [
               ] )

